#!/usr/bin/env python

from pylab import *
import MySQLdb as mdb
import os
import sys
import csv
import progressbar
from sets import Set

TIMEOUT = 900

###############
# TODO move this in different file at some point

def getFamily(path):
    tokens = path.split('/')
    return tokens[5]

class Job:
    def __init__(self, job_id, description, tlimit, mlimit, pb_set_id, args, binary, log):
        self.id = job_id
        self.description = description
        self.time_limit = tlimit
        self.memory_limit = mlimit
        self.problem_set_id = pb_set_id
        self.binary = binary
        self.log = log

    def __str__(self):
        return  "Job "+ str(self.id)+ ":\n" + \
         "    "+ self.description+ "\n" + \
         "     time_limit "+ str(self.time_limit)+ "+ mem_limit "+ str(self.memory_limit)+ "\n" +\
         "     problem_set "+ str(self.problem_set_id) + "\n" +\
         "     binary "+ self.binary+ "\n" +\
         "     log "+ self.log
        
    @staticmethod
    def getFromDB(cur, job_id):
        cur.execute("""SELECT
                id, description, time_limit, memory_limit,
                problem_set_id, args, timestamp, binary_path, log_path
                from Jobs where Jobs.id = %s;""", job_id)
        jobs = cur.fetchall()
        assert(len(jobs)==1)
        assert(jobs[0][0] == job_id)
        
        description = jobs[0][1]
        time_limit = int(jobs[0][2])
        memory_limit = int(jobs[0][3])
        problem_set_id = int(jobs[0][4])
        args = jobs[0][5]
        binary = jobs[0][7]
        log = jobs[0][8]
        return Job(job_id, description, time_limit,\
                   memory_limit, problem_set_id, args, binary, log)


class Result:
    Sat, Unsat, Unknown = "unsat", "sat", "unknown"

    @staticmethod
    def fromStr(string):
        if (string == "unsat"):
            return Result.Unsat
        if( string == "sat"):
            return Result.Sat;
        assert (string == "unknown")
        return Result.Unknown
    @staticmethod
    def str(res):
        if (res == Result.Sat):
            return "sat"
        if (res == Result.Unsat):
            return "unsat"
        if (res == Result.Unknown):
            return "unknown"

class Status:
     Timeout, Memout, Solved = "Timeout", "Memout", "Solved"

     @staticmethod
     def fromStr(string):
        if (string == "timeout"):
            return Status.Timeout
        if( string == "memout"):
            return Status.Memout
        assert (string == "solved")
        return Status.Solved

def write(txt):
    sys.stdout.write(txt)

def initializeProgressBar(max_value):
    """Set up the command-line progress bar with max_value
    """
    bar = progressbar.ProgressBar(maxval=max_value, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    return bar

def numLinesInFile(fname):
    """ Counts the number of lines in the given file.
    """
    with open(fname, 'rb') as f:
        for i, l in enumerate(f):
            pass
        return i + 1

class Problem:
    def __init__(self, pb_id, path):
        self.id = pb_id
        self.path = path

    @staticmethod
    def getFromDB(cur, problem_id):
        cur.execute("""SELECT id, path                 
                from Problems where Problems.id = %s;""", problem_id)
        pbs = cur.fetchall()
        assert(len(pbs)==1)
        assert(int(pbs[0][0]) == problem_id)
        return Problem(problem_id, pbs[0][1])

    def __hash__(self):
        return self.id
    def __str__(self):
        return "(id " + str(self.id) + ", path " + self.path + ")"
    def __eq__(self, other):
        if not isinstance(other, Problem):
            return False
        return self.id == other.id

class ProblemResult:
    def __init__(self, result, status, cpu_time, mem, exit_status):
        self.result = result
        self.status = status
        self.cpu_time = cpu_time
        self.memory = mem
        self.exit_status = exit_status

    def __str__(self):
        res = "(result: "+self.result + ", " \
              " status: "+self.status + ", " \
              " cpu_time: "+str(self.cpu_time) + ", " \
              " memory: "+str(self.memory) +")" \
              " exit_status: "+str(self.exit_status) + ")"
        return res
    def isSolved(self):
        return self.status == Status.Solved and \
               (self.result == Result.Sat or self.result == Result.Unsat)
    
    @staticmethod
    def getFromDB(cur, job_id, problem_id):
        cur.execute("""SELECT
                job_id, problem_id, result, cpu_time, memory, exit_status, status
                from JobResults where job_id = %s and problem_id = %s;""", (job_id, problem_id))
        jobs = cur.fetchall()
        assert(len(jobs)==1)
        assert(jobs[0][0] == job_id and jobs[0][1] == problem_id)
        
        result = Result.fromStr(jobs[0][2])
        cpu_time = float(jobs[0][3])
        memory = float(jobs[0][4])
        exit_status = int(jobs[0][5])
        status = Status.fromStr(jobs[0][6])
        
        return ProblemResult(result, status, cpu_time, memory, exit_status)

    @staticmethod
    def getFromDB(cur, job_id):
        cur.execute("""SELECT
                problem_id, result, cpu_time, memory, exit_status, status
                from JobResults where job_id = %s;""", (job_id))
        jresults = cur.fetchall()
        assert(len(jresults))

        pb_results = {}
        for jres in jresults:
            problem_id = int(jres[0])
            result = Result.fromStr(jres[1])
            cpu_time = float(jres[2])
            memory = float(jres[3])
            exit_status = int(jres[4])
            status = Status.fromStr(jres[5])
            pb_res = ProblemResult(result, status, cpu_time, memory, exit_status)
            pb_results[problem_id] = pb_res
        return pb_results
        

class ProblemSet:
    def __init__(self, pb_set_id, description, problems):
        self.id = pb_set_id
        self.description = description
        self.problems = problems

    @staticmethod
    def getFromDB(cur, problem_set_id):
        print "Getting problem set ", problem_set_id, " from data-base"
        cur.execute("""SELECT
                id, description from ProblemSets where id = %s;""", (problem_set_id))
        res = cur.fetchall()
        assert(len(res)==1)
        assert(res[0][0] == problem_set_id)
        description = res[0][1]

        cur.execute("""SELECT t1.problem_id, t2.path FROM ProblemSetToProblem t1
                       INNER JOIN Problems t2 ON t1.problem_id = t2.id
                       WHERE t1.problem_set_id = %s;""", (problem_set_id))
        res = cur.fetchall()
        assert (len(res))

        problems = []
        for r in res:
            pb_id = int(r[0])
            path = r[1]
            pb = Problem(pb_id, path)
            problems.append(pb)
        return ProblemSet(problem_set_id, description, problems)

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        if not isinstance(other, ProblemSet):
            return False
        return self.id == other.id
    
class CumulativeResult:
    def __init__(self):
        self.solved = 0
        self.cpu_time = 0.0

    def __str__(self):
        res = "(solved: "+self.solved + ", " \
              " cpu_time: "+self.cpu_time + ")"
        return res
        
    def addResult(self, res):
        if (res.status != Status.Solved):
            return
        if (res.result == Result.Unknown):
            return

        self.solved += 1
        self.cpu_time += res.cpu_time

class SolverToCumulativeResult:
    def __init__(self):
        self.solved = 0
        self.total_time = 0.0
        
    def addResult(self, result):
        if (result.isSolved):
            self.solved += 1
            self.total_time += result.cpu_time
            
class CummulativeResults:
    def __init__(self):
        self.table = {}
        
    def addSolverResult(self, solver, result):
        if not solver in self.table:
            self.table[solver] = CumulativeResult()
        self.table[solver].addResult(result)

    def getResult(self, solver):
        assert (solver in self.table)
        return self.table[solver]
    
    def getBest():
        max_solved = 0
        min_time = None
        best = None
        for sv in self.table:
            entry = self.table[sv]
            # first iteration
            if best == None:
                max_solved = entry.solved
                min_time = entry.total_time
                best = sv
                continue
            # clearly not better
            if entry.solved < max_solved:
                continue
            if entry.solved > max_solved:
                max_solved = entry.solved
                min_time = entry.total_time
                best = sv
                continue

            if entry.total_time < min_time:
                min_time = entry.min_time
                best = sv
        return best

def getBest(solverToRes, solvers):
    max_solved = 0
    min_time = None
    best = None

    for sv in solvers:
        assert (sv in solverToRes)
        entry = solverToRes[sv]
        # first iteration
        if best == None:
            max_solved = entry.solved
            min_time = entry.total_time
            best = sv
            continue
        # clearly not better
        if entry.solved < max_solved:
            continue
        if entry.solved > max_solved:
            max_solved = entry.solved
            min_time = entry.total_time
            best = sv
            continue
        
        if entry.total_time < min_time:
            min_time = entry.min_time
            best = sv
    return best
    

    
class SummaryTable:
    def __init__(self, families, problem_set):
        self.solvers = []
        self.table = {}
        for f in families:
            self.table[f] = {}
        self.fam_count = {}
        for f in families:
            self.fam_count[f] = 0
        for pb in problem_set.problems:
            fam = getFamily(pb.path)
            self.fam_count[fam] = self.fam_count[fam] + 1

    def addSolver(self, solver):
        if (not solver in self.solvers):
            self.solvers.append(solver)
    def addResult(self, sv, pb, res):
        fam = getFamily(pb.path)
        res2solv = self.table[fam]
        if not sv in res2solv:
            res2solv[sv] = CumulativeResult()
        res2solv[sv].addResult(res)

    def removeFamilies(self, fams):
        for fam in fams:
            del self.table[fam]
            del self.fam_count[fam]
        
    def csv(self):
        print "family, count, ",
        for solver in self.solvers:
            sv_name = solver
            print sv_name, " solved, ",
            print sv_name, " time, ",
        print ""
        
        for fam in self.table:
            print fam, ", ", self.fam_count[fam], ", ",
            sv2res = self.table[fam]
            for solver in self.solvers:
                result = sv2res[solver]
                print result.solved, ", ",
                print round(result.cpu_time,2), ", ",
            print ""

    def laTex(self, solvers):
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print "% Printing table for ",
        for sv in solvers:
            print sv, " ",
        print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        
        n = len(solvers)
        write("\\begin{tabular}{| l |")
        for i in range(n):
            write(" r r |")
        write("}\n")
        # headers
        write("\\cline{2-"+str(2*n+1)+"}\n")
        write("\\multicolumn{1}{c|}{}\n")
        for solver in solvers:
            write("& \\multicolumn{2}{c|}{\\textsf{"+solver+"}}\n")
        write ("\\\\\n")
        write("\\cline{2-"+str(2*n+1)+"}\n")
        write ("\\multicolumn{1}{l}{\\tiny\\textsf{set}}\n")
        for solver in solvers:
            write("& \\multicolumn{1}{r}{\\tiny\\textsf{solved}} & \\multicolumn{1}{r}{\\tiny\\textsf{time (s)}} \n")

        write ("\\\\\n\\hline\n")
        i= 0
        # total for each solver
        summary = SolverToCumulativeResult()
        
        for fam in sorted(self.table.keys()):
            if i % 2 == 0:
                write("\\rowcolor[gray]{.95}\n")
            i+=1
            # family
            fam_str = fam.replace("_","\_")
            write("\\textsf{"+fam_str+"}\n")
            # count
            write("("+str(self.fam_count[fam])+")\n")
            best = getBest(self.table[fam], solvers)
            # results
            sv_to_res = self.table[fam]
            for sv in solvers:
                res = sv_to_res.getResult(sv)
                if sv == best:
                    time = "\\textbf{"+str(round(res.cpu_time,2))+"} "
                    solved = "\\textbf{"+str(res.solved)+"} "
                else:
                    time = str(round(res.cpu_time,2))+" "
                    solved = str(res.solved)+" "
                write("& "+solved + " & " + time )
                summary.addSolverResult(sv, res)

            write("\n")
            write ("\\\\\n")
        write("\\hline\n")
        write("\\multicolumn{1}{c|}{}")
        # pick best sum
        best_total = summary.getBest()

        for sv in solvers:
            sum_res = summary.getResult(sv)
            if sv == best_total:
                time = "\\textbf{"+str(round(sum_res.total_time,2))+"} "
                solved = "\\textbf{"+str(sum_res.solved)+"} "
            else:
                time = str(round(sum_res.cpu_time,2))+" "
                solved = str(sum_res.solved)+" "
            write("& "+solved + " & " + time + "\n")
            
        write("\\\\\n\\cline{2-"+str(2*n+1)+"}\n")
        write("\\end{tabular}\n")
            
        
class SolverToResult:
    def __init__(self):
        self.table = {}

    def addResult(self, solver, result):
        assert (not solver in self.table)
        self.table[solver] = result

    def hasResult(self, solver):
        return (solver in self.table)
        
    def getResult(self, solver):
        assert (solver in self.table)
        return self.table[solver]

    def allAgree(self):
        result = None
        for sv in self.table:
            res = self.table[sv]
            if (res.status != Status.Solved):
                continue
            if (res.result == Result.Unknown):
                continue
            if (result == None):
                result = res.result
            if (result != res.result):
                return False
        return True

class ResultsTable:
    def __init__(self):
        self.table = {} # pb_id => SolverToResult
        self.solvers = Set()
        self.id_to_problem = {}
        self.problem_to_id = {}
        self.id_count = 0
        
    def newId(self):
        self.id_count += 1
        return self.id_count

    def getResult(self, solver, problem):
        assert solver in self.solvers
        sv_res = self.table[solver]
        pb_id = self.problem_to_id[problem]
        assert pb_id in sv_res
        return sv_res[pb_id]

    def addResult(self, problem, solver, result):
        if not (problem in self.table):
            self.table[problem] = SolverToResult()
        self.table[problem].addResult(solver, result)

    def addProblem(self, problem):
        if not problem in self.problem_to_id:
            pb_id = self.newId()
            self.problem_to_id[problem] = pb_id
            self.id_to_problem[pb_id] = problem

        return self.problem_to_id[problem]

    def addSolver(self, sv_config):
        if (not sv_config in self.solvers):
            self.solvers.add(sv_config)

    # adds results for the job to the results of the given solver
    def addJobToSolver(self, cur, solver, job_id):
        print "Adding job ", job_id, " to solver ", solver
        job = Job.getFromDB(cur,  job_id)
        problem_set = ProblemSet.getFromDB(cur, job.problem_set_id)

        print "Getting results for job ", job_id, "..."
        bar = initializeProgressBar(len(problem_set.problems))
        
        job_results = ProblemResult.getFromDB(cur, job_id)
        assert len(job_results) == len(problem_set.problems), "Missing problems from job"
        i = 0
        for pb in problem_set.problems:
            i = i + 1
            bar.update(i)
            pb_id = pb.id
            pb_local_id = self.addProblem(pb)
            pb_res = job_results[pb_id]
            self.addResult(pb_local_id, solver, pb_res)

    def addJobToSolverIncomplete(self, cur, solver, job_id):
        print "Adding job ", job_id, " to solver ", solver
        job = Job.getFromDB(cur,  job_id)
        problem_set = ProblemSet.getFromDB(cur, job.problem_set_id)

        print "Getting results for job ", job_id, "..."
        bar = initializeProgressBar(len(problem_set.problems))
        
        job_results = ProblemResult.getFromDB(cur, job_id)
        # assert len(job_results) == len(problem_set.problems), "Missing problems from job"
        i = 0
        for pb in problem_set.problems:
            i = i + 1
            bar.update(i)
            pb_id = pb.id
            pb_local_id = self.addProblem(pb)
            if pb_id in job_results:
                pb_res = job_results[pb_id]
                self.addResult(pb_local_id, solver, pb_res)
            
            
    def readRun(self, solver, bench_file, path):
        print "Reading run of solver ", solver, " from path ", path
        self.addSolver(solver)
        num_lines = numLinesInFile(bench_file)
        bench_count = 0
        bar = initializeProgressBar(num_lines)
        with open(bench_file, 'rb') as bfile:
            for bench in bfile:
                bench_count += 1
                bar.update(bench_count)
                bench = bench.strip()
                mangled_path = path + bench.replace("/","_")
                pb_id = self.addProblem(bench)
                res = ProblemResult.readResult(mangled_path)
                self.addResult(pb_id, solver, res)
                
            
    def printTable(self):
        print "problem, ",
        for solver in self.solvers:
            sv_name = solver
            print sv_name, "_result, ",
            print sv_name, "_cpu_time, ",
            print sv_name, "_memory, ",
            print sv_name, "_status, ",
            print sv_name, "_exit_status,",
        print ""
        
        for pb in self.id_to_problem:
            assert pb in self.table
            path = self.id_to_problem[pb].path
            restosolv = self.table[pb]
            print path, ", ",
            for solver in self.solvers:
                result = restosolv.getResult(solver)
                print result.result, ", ",
                print result.cpu_time, ", ",
                print result.memory, ", ",
                print result.status, ", ",
                print result.exit_status, ", ",
            print ""

    def uniqueSolves(self, sv1, sv2):
        errors = []
        sv1_solves = []
        sv2_solves = []
        for pb_id in self.table:
            problem = self.id_to_problem[pb_id]
            sv2res = self.table[pb_id]
            res1 = sv2res.getResult(sv1)
            res2 = sv2res.getResult(sv2)
            if (res1.result == res2.result):
                continue
            if (res1.result == Result.Sat or
                res1.result == Result.Unsat):
                assert (res2.result == Result.Unknown)
                sv1_solves.append((problem, res1))
            else:
                assert (res1.result == Result.Unknown)
                assert (res2.result == Result.Sat or
                        res2.result == Result.Unsat)
                sv2_solves.append((problem, res2))
                    
        print "Unique solves: ", sv1, " ", len(sv1_solves)
        for p in sv1_solves:
            print p[0].path, " ", p[1].cpu_time, " ", p[1].result

        print "Unique solves: ", sv2, " ", len(sv2_solves)
        for p in sv2_solves:
            print p[0].path, " ", p[1].cpu_time, " ", p[1].result

    def uniqueSolvesIncomplete(self, sv1, sv2):
        errors = []
        sv1_solves = []
        sv2_solves = []
        for pb_id in self.table:
            problem = self.id_to_problem[pb_id]
            sv2res = self.table[pb_id]
            if not sv2res.hasResult(sv1) or not sv2res.hasResult(sv2):
                continue
            
            res1 = sv2res.getResult(sv1)
            res2 = sv2res.getResult(sv2)
            if (res1.result == res2.result):
                continue
            if (res1.result == Result.Sat or
                res1.result == Result.Unsat):
                assert (res2.result == Result.Unknown)
                sv1_solves.append((problem, res1))
            else:
                assert (res1.result == Result.Unknown)
                assert (res2.result == Result.Sat or
                        res2.result == Result.Unsat)
                sv2_solves.append((problem, res2))
                    
        print "Unique solves: ", sv1, " ", len(sv1_solves)
        for p in sv1_solves:
            print p[0].path, " ", p[1].cpu_time, " ", p[1].result

        print "Unique solves: ", sv2, " ", len(sv2_solves)
        for p in sv2_solves:
            print p[0].path, " ", p[1].cpu_time, " ", p[1].result
            

    def checkCorrectness(self):
        errors = []
        for pb_id in self.table:
            problem = self.id_to_problem[pb_id]
            sv2res = self.table[pb_id]
            result = None
            for sv in sv2res.table:
                res = sv2res.getResult(sv)
                if (res.status != Status.Solved):
                    continue
                if (res.result == Result.Unknown):
                    errors.append((problem, sv))
                    continue
                    
                if (res.result == Result.Sat or res.result == Result.Unsat):
                    if (result == None):
                        result = res.result
                        
                    if (result != res.result):
                        print "Different answer on problem: "
                        print problem
                        for s in sv2res.table:
                            print s, ": ", sv2res.getResult(s).result, " time: ", sv2res.getResult(s).cpu_time
                    
        print "\nPotential crashes: "
        for er in errors:
            print "    ", er[1], " on ", er[0]

    def removeFamilies(self, fams):
        remove = Set()
        for pb_id in self.table:
            problem = self.id_to_problem[pb_id]
            for fam in fams:
                if (problem.path.find(fam) != -1):
                    remove.add((pb_id, problem))
                    break
                
        for (pb_id, pb) in remove:
            del self.table[pb_id]
            del self.id_to_problem[pb_id]
            del self.problem_to_id[pb]
            
    def removeCrashes(self):
        errors = []
        for pb_id in self.table:
            problem = self.id_to_problem[pb_id]
            sv2res = self.table[pb_id]
            crash = False
            for sv in sv2res.table:
                res = sv2res.getResult(sv)
                if (res.result == Result.Unknown and
                    res.status == Status.Solved and
                    res.cpu_time < 290):
                    crash = True
            if crash:
                errors.append((pb_id, problem))
        print "\nPotential crashes: ",len(errors)
        for er in errors:
            print "    ", er[1], " on ", er[0]
        print "\nRemoving crashes.."
        for (pb_id, pb) in errors:
            del self.table[pb_id]
            del self.id_to_problem[pb_id]
            del self.problem_to_id[pb]
            
    def getCactusData(self, solvers):
        data = {}
        for sv in solvers:
            data[sv] = []
        for pb in self.table:
            sv2res = self.table[pb]
            for sv in solvers:
                sv_res = sv2res.table[sv]
                if (sv_res.isSolved()):
                    data[sv].append(float(sv_res.cpu_time))
        for sv in solvers:
            data[sv] = sort(data[sv])
        return data
                
            
    def csv(self):
        print "problem, ",
        for solver in self.solvers:
            sv_name = solver
            print sv_name, " result, ",
            print sv_name, " cpu_time, ",
            print sv_name, " mem, ",
            print sv_name, " status, ",
            print sv_name, " exit, ", 
        print ""
        
        for pb in self.id_to_problem:
            assert pb in self.table
            path = self.id_to_problem[pb]
            restosolv = self.table[pb]
            print path, ", ",
            for solver in self.solvers:
                result = restosolv.getResult(solver)
                print result.result, ", ",
                print result.cpu_time, ", ",
                print result.mem, ", ",
                print result.status, ", ",
                print result.exit_status, ", ",
            print ""
            
    def fromCSV(self, path):
        fields = ["result", "cpu_time", "mem", "status", "exit"]
        num_lines = numLinesInFile(path)
        bar = initializeProgressBar(num_lines - 1)
        num_problems = 0
        solvers = Set ()
        with open(path, 'rb') as csvfile:
            reader = csv.DictReader(csvfile)
            headers = reader.next()
            for h in headers:
                tokens = h.split()
                if (len(tokens) == 0):
                    continue
                if (len(tokens) == 1):
                    assert (tokens[0] == "problem")
                    continue
                solver = tokens[0]
                assert (len(tokens) == 2 and tokens[1] in fields)
                self.addSolver(solver)
                solvers.add(solver)
                
            for row in reader:
                num_problems +=1
                bar.update(num_problems)
                problem = row["problem"]
                pb_id = self.addProblem(problem)
                for sv in solvers:
                    result = row["  "+sv+"  "+"result"].strip()
                    cpu_time = float(row["  "+sv+"  "+"cpu_time"])
                    mem = float(row["  "+sv+"  "+"mem"])
                    status = row["  "+sv+"  "+"status"].strip()
                    exit_status = int(row["  "+sv+"  "+"exit"])
                    res = ProblemResult(result, status, cpu_time, 0, mem, exit_status)
                    self.addResult(pb_id, sv, res)

    def mkSummary(self, families, problem_set):
        sum_table = SummaryTable(families, problem_set)

        for  sv in self.solvers:
            sum_table.addSolver(sv)
        for pb_id in self.table:
            pb = self.id_to_problem[pb_id]
            sv2res = self.table[pb_id]
            for sv in sv2res.table:
                res = sv2res.getResult(sv)
                sum_table.addResult(sv, pb, res)
        return sum_table
            
            
    def summary(self):
        crashes = {}
        cres = {}
        errors = []

        for solver in self.solvers:
            cres[solver] = CumulativeResult()
            crashes[solver] = []

            
        for pb in self.id_to_problem:
            assert pb in self.table
            path = self.id_to_problem[pb]
            restosolv = self.table[pb]
            # skipping problems the solvers disagree on
            if (not restosolv.allAgree()):
                errors.append((path, restosolv))
                continue

            for sv in self.solvers:
                result = restosolv.getResult(sv)
                cres[sv].addResult(result)
                if (result.status == Status.Solved and result.result == "Unknown"):
                    # probably a crash
                    crashes[sv].append(path)

        for sv in self.solvers:
            print sv, " solved, ",
            print sv, " time, ",
            print sv, " crashes, ",
        print ""
        
        for sv in self.solvers:
            print cres[sv].solved, ", ",
            print cres[sv].cpu_time, ", ",
            print len(crashes[sv]), ", ",
        print ""

        print "Solvers do not agree on the following problems "
        for (pb, sv2res) in errors:
            print "On problem ", pb 
            for sv in self.solvers:
                print sv, "   ", sv2res.getResult(sv).result
            print ""

    def getScatterResults(self):
        results = []
        for pb in self.id_to_problem:
            solver2res = self.table[pb]
            result = Result.Unknown
            assert (solver2res.allAgree())
            for sv in solver2res.table:
                res = solver2res.getResult(sv)
                if (res.result == Result.Unknown):
                    continue
                if (result == Result.Unknown):
                    result = res.result
                else:
                    continue
                    assert (result == res.result)
            results.append(result)
        return results

    def getScatterResultsIncomplete(self, sv1, sv2):
        results = []
        for pb in self.id_to_problem:
            if not pb in self.table:
                continue;
            solver2res = self.table[pb]
            if not solver2res.hasResult(sv1) or not solver2res.hasResult(sv2):
                continue
            result = Result.Unknown
            assert (solver2res.allAgree())
            for sv in solver2res.table:
                res = solver2res.getResult(sv)
                if (res.result == Result.Unknown):
                    continue
                if (result == Result.Unknown):
                    result = res.result
                else:
                    continue
                    assert (result == res.result)
            results.append(result)
        return results
    
    
    def getScatterData(self, xsolver, ysolver, timeout):
        xs = []
        ys = []
        for pb in self.id_to_problem:
            solver2res = self.table[pb]
            x_res = solver2res.getResult(xsolver)
            y_res = solver2res.getResult(ysolver)
            
            x_point = timeout
            y_point = timeout
            if (not solver2res.allAgree()):
                print "Solvers do not agree on ",pb, ": ", self.id_to_problem[pb]
                for sv in solver2res.table:
                    print "    ",sv,": ",solver2res.table[sv] 
                assert(false)
            if (x_res.status == Status.Solved and x_res.result != "Unknown"):
                x_point = x_res.cpu_time

            if (y_res.status == Status.Solved and y_res.result != "Unknown"):                
                y_point = y_res.cpu_time

            xs.append(x_point)
            ys.append(y_point)
        return (xs, ys)

    def getScatterDataIncomplete(self, xsolver, ysolver, timeout):
        xs = []
        ys = []
        for pb in self.id_to_problem:
            if not pb in self.table:
                continue
            solver2res = self.table[pb]
            if not solver2res.hasResult(xsolver) or not solver2res.hasResult(ysolver):
                continue
            x_res = solver2res.getResult(xsolver)
            y_res = solver2res.getResult(ysolver)
            
            x_point = timeout
            y_point = timeout
            if (not solver2res.allAgree()):
                print "Solvers do not agree on ",pb, ": ", self.id_to_problem[pb]
                for sv in solver2res.table:
                    print "    ",sv,": ",solver2res.table[sv] 
                assert(false)
            if (x_res.status == Status.Solved and x_res.result != "Unknown"):
                x_point = x_res.cpu_time

            if (y_res.status == Status.Solved and y_res.result != "Unknown"):                
                y_point = y_res.cpu_time

            xs.append(x_point)
            ys.append(y_point)
        return (xs, ys)

    
def makeScatterPlot(title, x_data, y_data, results, x_label, y_label, limit):
    font = {'family' : 'sans-serif',
            'weight' : 'normal',
            'size'   : 22}
    rc('font', **font)
    figure(figsize=(7, 7), dpi=80)
    #ax = fig.add_subplot(1,1,1) # one row, one column, first plot
    ax = gca()
    gcf().subplots_adjust(bottom=0.18)
    gcf().subplots_adjust(left=0.18)

    for k, spine in ax.spines.items():  #ax.spines is a dictionary
        spine.set_zorder(-1)
        spine.set_position(('outward', 10))
        
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')

    if limit == None:
        limit = max(max(x_data), max(y_data))
        low_limit = min(min(x_data), min(y_data))
    else:
        low_limit = 0
        
    xlim(low_limit, limit)
    ylim(low_limit, limit)
    
    x_diagonal = range(int(limit))
    plot(x_diagonal, x_diagonal, color="grey", linestyle="-", linewidth="0.3")
    grid(True)
    # ax.set_tick_params('both','minor','off')

    if (results == None):
        scatter(x_data, y_data, color="red", marker="^", clip_on=False)
    else:
        # split Failure and Success
        xsuccess, xfailed, xunknown = [], [], []
        ysuccess, yfailed, yunknown = [], [], []

        for i, res in enumerate(results):
            if (res == Result.Unsat):
                xfailed.append(x_data[i])
                yfailed.append(y_data[i])
            elif res == Result.Sat:
                xsuccess.append(x_data[i])
                ysuccess.append(y_data[i])
            else:
                xunknown.append(x_data[i])
                yunknown.append(y_data[i])
        scatter(xfailed, yfailed, color="red", marker="x", clip_on=False)
        scatter(xsuccess, ysuccess, color="green", marker="^", clip_on=False)
        scatter(xunknown, yunknown, color="black", marker="o", clip_on=False)
    
    # ax.set_title(solverx+" vs "+solvery)
    
    xlabel(x_label, fontsize=22, fontweight="bold")
    ylabel(y_label, fontsize=22, fontweight="bold")
    
    savefig(title+".png")
    print "Created "+title+".png"
    xscale('symlog')
    yscale('symlog')
    
    savefig(title+"_log.png")
    print "Created "+title+"_log.png"

def makeCactusPlot(solvers,title, xlimit):
    font = {'family' : 'sans-serif',
            'weight' : 'normal',
            'size'   : 15}
    rc('font', **font)
    figure(figsize=(10, 7), dpi=80)
    #ax = fig.add_subplot(1,1,1) # one row, one column, first plot
    ax = gca()
    gcf().subplots_adjust(bottom=0.18)
    gcf().subplots_adjust(left=0.22)


    for k, spine in ax.spines.items():  #ax.spines is a dictionary
        spine.set_zorder(-1)
        spine.set_position(('outward', 10))
        
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')

    grid(True)
    
    xmax_limit = 0
    ymax_limit = 300
    for sv in solvers:
        sv_data = solvers[sv]
        xmax_limit = max(xmax_limit, len(sv_data))
        ymax_limit = max(ymax_limit, sv_data[-1])

    xlim(xlimit, xmax_limit)
    ylim(0, ymax_limit)

    
    colors = ['r', 'b', 'y', 'm', 'c','g','k','brown','orange']
    markers = ['v', 'o', 's', '^', '*','x','v','.','s']
    assert (len(colors) >= len(solvers))
    
    scatters = []
    names = []
    i = 0
    for sv in solvers:
        sv_data = solvers[sv]
        num_solved = range(0,len(sv_data))
        if xlimit != None:
            sv_data = sv_data[xlimit:]
            num_solved = num_solved[xlimit:]
        sc = scatter(num_solved, sv_data, color=colors[i], s=12, marker=markers[i], clip_on=False)
        plot(num_solved, sv_data, color=colors[i])
        scatters.append(sc)
        names.append(sv)
        i+=1
        
    legend(tuple(scatters),
           tuple(names),
           scatterpoints=1,
           loc='upper left',
           ncol=1)

    xlabel("# solved", fontsize=22, fontweight="bold")
    ylabel("time (s)", fontsize=22, fontweight="bold")
    
    savefig(title+".png")
    print "Created plot ",title+".png"
    # xscale('symlog')
    # yscale('symlog')
    
    savefig(title+"_log.png")
    print "Created plot ",title+"_log.png"

def makeBenchmarkBar(data, results, timeout, title):
    bins = [10*x for x in range(31)]
    # bins = [1, 2, 3, 4100,125,150,160,170,180,190,200,210,220,230,240,250,275,300]
    # the histogram of the data with histtype='step'
    data_fail = []
    data_suc = []
    data_unk = []
    for i, res in enumerate(results):
        if (res == Result.Unsat):
            data_fail.append(data[i])
        elif (res == Result.Sat):
            data_suc.append(data[i])
        else:
            data_unk.append(data[i])
            
    n, bins, patches = hist([data_fail, data_suc, data_unk],bins, \
                                label = ["unsat", "sat", "unknown"], \
                                color = ['red', 'green', 'gray'], histtype='barstacked', rwidth=0.4)
    
    xlabel("runtime", fontsize=14, fontweight="bold")
    ylabel("number benchmarks", fontsize=14, fontweight="bold")
    legend(title=title)
    #
    # now we create a cumulative histogram of the data
    #
    grid(True)

    print "Creating figure ", title, "_histogram.png"
    savefig(title+"_histogram.png")
    figure()

