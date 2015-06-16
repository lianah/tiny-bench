#!/usr/bin/env python

import MySQLdb as mdb
import sys
import datetime
import argparse
import os
import os.path
import subprocess
import benchmarkingUtils as bu

parser = argparse.ArgumentParser(description='Adds a job the the database')
parser.add_argument('-b', '--binary',type=str,
                    help='path to binary', required=True)
parser.add_argument('-p', '--problemset', type=int,
                    help='problem set id', required=True)
parser.add_argument('-t', '--time', type=float,
                    help='maximum time per benchmark (s)', default=100.0)
parser.add_argument('-m', '--memory', type=float,
                    help='maximum time per benchmark (MB)', default=1500.0)
parser.add_argument('-d', '--description', type=str,
                    help='A longer description of the job', default="")
parser.add_argument('--manual-start', help='don\'t start the runner automatically',
                   action="store_true")
parser.add_argument('-l', '--logpath', type=str, required=True,
                     help='path where to store log files',
                    default=[])
parser.add_argument('-a', '--arguments', nargs=argparse.REMAINDER,
                     help='command line arguments for binary',
                    default=[])

args = parser.parse_args()

binary_path = args.binary
#assert os.path.isfile(binary_path), "Cannot find file:" + binary_path

problem_set = args.problemset
arguments = ' '.join(args.arguments)
time_limit = args.time
memory_limit = args.memory
description = args.description
log_path = args.logpath
auto_start = not(args.manual_start)

print "auto_start: ", auto_start


def addJob(cur, description, time_limit, memory_limit,
           problem_set, arguments, binary_path, log_path):
    cur.execute("""INSERT into Jobs
                (id, description, time_limit, memory_limit,
                 problem_set_id, args, timestamp, binary_path, log_path)
                VALUES(default, %s, %s, %s, %s,
                 %s, default, %s, %s);""",
                (description, time_limit, memory_limit,
                 problem_set, arguments, binary_path, log_path))
    numrows = int(cur.rowcount)
    assert numrows == 1
    cur.execute("""SELECT LAST_INSERT_ID();""")
    job_id_res = cur.fetchone()
    assert len(job_id_res) == 1
    return job_id_res[0]

def addToQueue(cur, job_id, problem_set_id):
    cur.execute("""insert into Queue (job_id, problem_id)
                   select %s, Problems.id
                   FROM Problems INNER JOIN ProblemSetToProblem
                   ON Problems.id=ProblemSetToProblem.problem_id
                   WHERE ProblemSetToProblem.problem_set_id=%s;""",
                (job_id, problem_set_id))
    return int(cur.rowcount)

(server, user, password, database) = bu.loadConfiguration()
con = mdb.connect(server, user, password, database);
with con:
    cur = con.cursor()
    job_id = addJob(cur, description, time_limit, memory_limit,
                    problem_set, arguments, binary_path, log_path)

    print "Adding job:",
    bu.printJobWithId(cur, job_id)
    added = addToQueue(cur, job_id, problem_set)
    print "Inserted Queue ", added, " problem into Queue"
con.close()


# # starting the runner
# if auto_start:
#     print "Starting runner for job", job_id, "..."
#     subprocess.call(["./runner.py", str(job_id)])
# else:
#     print "Start a runner with job ID", job_id
#     subprocess.call(["./createpbs.py", str(job_id), "jobtorun.pbs"])
