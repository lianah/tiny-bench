#!/usr/bin/env python

import MySQLdb as mdb
import sys
import os
import argparse
import subprocess
import re
import benchmarkingUtils as bu
import platform, hashlib
import paramiko

RUN_SOLVER="runsolver"

(server, user, password, database) = bu.loadConfiguration()

parser = argparse.ArgumentParser(description='Run a job.')
parser.add_argument('JobID', type=int, help='id of the job you want to work on')
parser.add_argument('-v', '--verbosity', type=int,
                    help='how verbose the script should be', default=0)
parser.add_argument('-m', '--machine', type=str,
                    help='machine to run process remotely on', default="dkr10")

args = parser.parse_args()

job_id = args.JobID
verbosity = args.verbosity
machine = args.machine

# PID : combination of both system name and current process ID
pid = os.getpid()
pid += (int(hashlib.md5(platform.node()).hexdigest(), 16) % 10000) * 100000


def storeProblemResult(problem_id, run_time, memory, result, exit_status, status):
    problem_result_id=None
    con = mdb.connect(server, user, password, database);
    with con:
        cur = con.cursor()
        cur.execute("""insert into JobResults (job_id, problem_id, cpu_time, memory, result, exit_status, status, machine)
                    VALUES(%s, %s, %s, %s, %s, %s, %s, %s);""",
                    (job_id, problem_id, run_time, memory, result, exit_status, status, machine))
        cur.execute("""select id from JobResults
                       where job_id=%s and problem_id=%s;""",
                    (job_id, problem_id));
        probResults = cur.fetchall()
        assert len(probResults) == 1, "Failed to store problem result: Problem may already be in the database"
        assert len(probResults[0])== 1
        problem_result_id = probResults[0][0]

        #remove from Queue
        cur.execute("""DELETE FROM Queue
                       WHERE problem_id=%s AND job_id=%s;""",
                    (problem_id, job_id))
        numrows = int(cur.rowcount)
        assert numrows == 1, "Failed to delete job from queue: Concurrency error?"
    con.close()
    assert problem_result_id != None
    return problem_result_id

def storeError(job_result_id):
    con = mdb.connect(server, user, password, database);
    with con:
        cur = con.cursor()
        cur.execute("""update JobResults set result='error' where id='%s'""",
                    (job_result_id));
        print "Setting the result of ", job_result_id, "to 'error'"
        print "touched", int(cur.rowcount)
    con.close()


def runRemoteProcess(problem_id, problem_path, err_log, out_log, runsolver_log):
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(
        paramiko.AutoAddPolicy())
    ssh.connect(machine, username='lianah')

    run_cmd = RUN_SOLVER + " -C " + str(time_limit) + \
        " -S " + str(mem_limit) + " -w " + runsolver_log + \
        " " + " -o " + out_log + \
        " " + binary_path + " " + args + " " + problem_path

    if verbosity > 0:
        print "    " + run_cmd

    job_log_path = log_path + str(job_id)
    check_dir_cmd = 'if test ! -d ' + job_log_path + ' ; then mkdir ' + job_log_path + '; fi'
    stdin, stdout,stderr = ssh.exec_command(check_dir_cmd)
    out = stdout.readlines()
    err = stderr.readlines()
    if (len(out) or len(err)):
        print "ERROR running ", check_dir_cmd
        print out
        print err
        exit

    stdin, stdout,stderr = ssh.exec_command("cd ~/benchmarking/simpler-benchmarking; " + run_cmd)
    err = stderr.readlines()
    
    if (len(err)):
        print "ERROR running ", run_cmd
        for e in err:
            print e
            
    out = stdout.readlines()
    for l in out:
        print l

    cmd = 'echo "[Benchmark Path] ' + problem_path + '" >> ' + out_log
    stdin, stdout,stderr = ssh.exec_command(cmd);

    out = stdout.readlines()
    err = stderr.readlines()
    assert (not len(out) and not len(err))

    cmd = '/users/lianah/benchmarking/simpler-benchmarking/grab-results.py -o ' + out_log + ' -w ' + runsolver_log
    stdin, stdout,stderr = ssh.exec_command(cmd);
    assert (not len(stderr.readlines()))

    out = stdout.readlines()
    (cpu, memory, result, exit_status, status) = (None, None, None, None, None)
    for line in out:
        if (line.startswith("memory")):
            memory = float (line [len("memory") + 1:].strip())
        elif (line.startswith("cpu")):
            cpu = float (line [len("cpu") + 1:].strip())
        elif (line.startswith("exit_status")):
            exit_status = int (line [len("exit_status") + 1:].strip())
        elif (line.startswith("result")):
            result = line[len("result")+1:].strip()
        elif (line.startswith("status")):
            status = line[len("status")+1:].strip()

    assert (cpu != None and \
            memory != None and\
            result != None and\
            exit_status != None and\
            status != None)

    return (cpu, memory, result, exit_status, status)
            

def runProblem(p):
    (problem_id, problem_path) = p
    print "Running", problem_path, "..."

    # error output log
    err_log=bu.genLogPath(log_path, job_id, problem_id, ".err")
    out_log=bu.genLogPath(log_path, job_id, problem_id, ".out")
    runlim_log=bu.genLogPath(log_path, job_id, problem_id, ".w")

    (run_time, memory, result, exit_status, status) = runRemoteProcess(problem_id, problem_path, err_log, out_log, runlim_log)

    print "Result ", result, "(", status, ")"
    print "in (", run_time, "seconds,", memory, " MB) with exit status ", exit_status
    
    job_result_id = storeProblemResult(problem_id, run_time, memory, result, exit_status, status)
    return job_result_id

def popQueue():
    con = mdb.connect(server, user, password, database);
    emp, problem = None, None
    with con:
        cur = con.cursor()
        cur.execute("""select count(*) from Queue
                       where runner_pid IS NULL and job_id=%s;""",
                    (job_id))
        res = cur.fetchall()
        assert len(res) == 1
        assert len(res[0]) == 1
        qs = res[0][0]
        if qs <= 0:
            emp = True
        else:
            cur.execute("""update Queue set runner_pid=%s
                           where Queue.runner_pid IS NULL and Queue.job_id=%s
                           limit 1;""",
                        (pid, job_id))
            cur.execute("""select Problems.id, Problems.path
                           FROM Problems INNER JOIN Queue ON
                           Problems.id=Queue.problem_id
                           WHERE Queue.job_id=%s and Queue.runner_pid=%s
                           limit 1;""",
                        (job_id, pid))
            problems = cur.fetchall()
            if len(problems) >= 1:
                problem = problems[0]
    con.close()
    return (emp, problem)

# Grab globals for job
print "Running job", job_id, "with process id", pid
(time_limit, mem_limit, binary_path, args, log_path) = (None, None, None, None, None)
con = mdb.connect(server, user, password, database);
with con:
    cur = con.cursor()
    cur.execute("""SELECT
                   time_limit, memory_limit, binary_path, args, log_path
                   FROM Jobs WHERE Jobs.id=%s;""",
                (job_id))
    res = cur.fetchall()
    assert len(res) == 1

    (time_limit, mem_limit, binary_path, args, log_path) = res[0]
    time_limit=int(time_limit)
    mem_limit=int(mem_limit)
con.close()

print "With time_limit=",time_limit, ", mem_limit=",mem_limit,
print ", args=",args,", binary_path=",binary_path,
print ", log_path=",log_path


#main loop
while True:
    (emp, problem) = popQueue()
    if emp == True:
        print "Done"
        break
    elif problem != None:
        runProblem(problem)
