import sys
import os

USER_FILE="./config/user"
PASSWORD_FILE="./config/password"
HOST_FILE="./config/host"
DATABASE_FILE="./config/database"

def echoFile(fileName):
    try:
        f = open(fileName, 'r');
        contents = f.read()
        return contents.strip()
    except IOError:
        sys.exit("Could not open " + fileName+". Make sure this exists and is readable.")

def loadConfiguration():
    server=echoFile(HOST_FILE)
    user=echoFile(USER_FILE)
    password=echoFile(PASSWORD_FILE)
    database=echoFile(DATABASE_FILE)
    return (server, user, password, database)

def printJobWithId(cur, job_id):
    printJob(selectJob(cur, job_id))

def selectJob(cur, job_id):
    cur.execute("""SELECT
                id, description, time_limit, memory_limit,
                problem_set_id, args, timestamp, binary_path, log_path
                from Jobs where Jobs.id = %s;""", job_id)
    jobs = cur.fetchall()
    assert(len(jobs)==1)
    assert(jobs[0][0] == job_id)
    return jobs[0]


def printJob(j):
    job_id = j[0]
    description = j[1]
    time_limit = j[2]
    memory_limit = j[3]
    problem_set_id = j[4]
    arguments = j[5]
    timestamp = j[6]
    binary_path = j[7]
    log_path = j[8]

    print job_id, ":", description
    print " ", "Problem set id: ", problem_set_id
    print " ", timestamp.isoformat(), time_limit, memory_limit
    print " ", binary_path, arguments
    print " log: ", log_path

def createPath(path):
    if not os.path.isdir(path):
        os.mkdir(path)

def genLogPath(log_path, job_id, problem_id, ext):
    return log_path+str(job_id)+"/"+str(problem_id)+ext
        

