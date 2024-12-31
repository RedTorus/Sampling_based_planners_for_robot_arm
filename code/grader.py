import numpy as np
import pdb
import argparse
import subprocess # For executing c++ executable
import pandas as pd
from timeit import default_timer as timer

plannerList = ["RRT", "RRTCONNECT", "RRTSTAR", "PRM"]

###############################################################
################### Util Functions Below ######################

def convertPIs(aString):
    """ Input: A comma seperated string like "pi/2,pi/4,pi/2,pi/4,pi/2,"
            or 1.2,4,5.3 etc
    Output: string replacing the pis 1.57079,...,...
    """
    if aString[-1] == ",":  # Remove training comma if there is one
        aString = aString[:-1]
    aString = aString.replace("pi", "3.141592") # Replace pi with 3.14... if needed
    vecOfStrings = aString.split(",")
    ans = []
    for anExpression in vecOfStrings:
        ans.append(str(eval(anExpression))) # Evaluate expressions if needed
    return ans

###############################################################
################### Main Functions Below ######################


# ["map1.txt", "1.570796,0.785398,1.570796,0.785398,1.570796",
#                                 # "0.392699,2.356194,3.141592,2.8274328,4.712388"],
#                                 "0.392699,2.356194,3.141592,2.8274328,4.712388"],
#             ["map2.txt", "0.392699,2.356194,3.141592",
#                                 "1.570796,0.785398,1.570796"],

def graderMain(executablePath, gradingCSV):

    # problems = [["map1.txt", "1.570796,0.785398,1.570796,0.785398,1.570796",
    #                             # "0.392699,2.356194,3.141592,2.8274328,4.712388"],
    #                             "0.392699,2.356194,3.141592,2.8274328,4.712388"],
    #         ["map2.txt", "0.392699,2.356194,3.141592",
    #                             "1.570796,0.785398,1.570796"]]

    # problems = [
    #         ["map2.txt", "1.738101,0.392955,0.478313,3.695973", "0.205897,2.476178,2.082369,0.166658"],
    #         ["map2.txt", "0.923738,1.099271,1.228511,5.380957", "0.792536,2.722955,5.032667,1.163456"],
    #         ["map2.txt", "1.496942,1.928107,3.119629,5.617534", "0.624424,2.074921,3.403268,0.292450"],
    #         ["map2.txt", "0.901069,2.187607,4.691609,6.088633", "1.821428,6.022235,2.830105,3.422870"],
    #         ["map2.txt", "0.686219,2.829659,2.962882,4.023547", "0.605440,2.605078,0.580391,4.301691"],
    #         ["map2.txt", "0.711608,3.389074,0.345863,5.989326", "0.960064,2.680343,2.311031,0.985693"],
    #         ["map2.txt", "1.080360,2.221731,2.865437,1.612382", "0.002021,0.363770,1.091501,4.316685"],
    #         ["map2.txt", "0.502838,1.659344,0.174526,5.158189", "0.474594,0.530737,2.135372,0.452818"],
    #         ["map2.txt", "1.047488,1.219040,3.234664,6.140207", "0.343790,2.641898,6.070597,1.711953"],
    #         ["map2.txt", "1.630170,2.781529,1.060634,4.591730", "1.251509,2.102154,3.065524,2.798367"],
    #         ["map2.txt", "0.591074,1.423401,2.456847,3.080740", "0.962081,2.635202,5.083832,1.430397"],
    #         ["map2.txt", "0.905398,2.751766,6.208815,1.261215", "1.316464,2.843499,3.723832,1.162286"],
    #         ["map2.txt", "0.835794,1.545433,3.829067,0.805480", "1.266769,0.922631,4.807138,5.516227"],
    #         ["map2.txt", "1.391913,1.218437,3.380623,3.336154", "0.858409,3.648787,0.853047,3.952140"],
    #         ["map2.txt", "1.643638,2.320623,2.684987,1.341215", "0.113686,1.150052,0.176047,2.856480"],
    #         ["map2.txt", "1.674122,0.596322,3.105708,0.058633", "0.496465,6.114873,3.267507,0.423600"],
    #         ["map2.txt", "1.290643,2.044831,0.571514,4.217606", "0.517610,2.391783,2.535478,1.014074"],
    #         ["map2.txt", "1.156436,1.173089,3.640363,0.471278", "1.377275,2.784672,6.209863,5.961038"],
    #         ["map2.txt", "1.840387,6.175411,5.348832,0.700952", "0.925034,2.910665,0.239961,3.951177"],
    #         ["map2.txt", "1.865864,1.719841,6.254910,0.774910", "1.200433,1.961238,0.603995,5.688892"],
    #         ]
    
    problems = [["map2.txt", "1.738101,0.392955,0.478313,3.695973", "0.205897,2.476178,2.082369,0.166658"],
            ["map2.txt", "1.738101,0.392955,0.478313,3.695973", "0.205897,2.476178,2.082369,0.166658"],
            ["map2.txt", "1.738101,0.392955,0.478313,3.695973", "0.205897,2.476178,2.082369,0.166658"],
            ["map2.txt", "1.738101,0.392955,0.478313,3.695973", "0.205897,2.476178,2.082369,0.166658"],
            ["map2.txt", "1.738101,0.392955,0.478313,3.695973", "0.205897,2.476178,2.082369,0.166658"],
            ["map2.txt", "1.738101,0.392955,0.478313,3.695973", "0.205897,2.476178,2.082369,0.166658"],
            ["map2.txt", "1.738101,0.392955,0.478313,3.695973", "0.205897,2.476178,2.082369,0.166658"],
            ["map2.txt", "1.738101,0.392955,0.478313,3.695973", "0.205897,2.476178,2.082369,0.166658"],
            ["map2.txt", "1.738101,0.392955,0.478313,3.695973", "0.205897,2.476178,2.082369,0.166658"],
            ["map2.txt", "1.738101,0.392955,0.478313,3.695973", "0.205897,2.476178,2.082369,0.166658"],]


    # problems =[["map2.txt", "0.138026,2.457708,1.506632", "0.986330,1.934048,3.974117"],
    #         ["map2.txt", "1.486140,1.339818,4.568666", "1.696222,0.306771,3.167160"],
    #         ["map2.txt", "0.726420,1.431127,1.275730", "0.978253,1.442233,3.490892"],
    #         ["map2.txt", "0.319935,2.168653,4.922019", "1.595666,0.586004,2.323059"],
    #         ["map2.txt", "1.547585,5.190207,0.117080", "1.677945,2.133129,1.384947"],
    #         ["map2.txt", "1.705660,1.167371,0.708549", "0.839724,2.533624,6.203962"],
    #         ["map2.txt", "0.114671,2.153189,2.257266", "1.497620,0.495936,1.507648"],
    #         ["map2.txt", "1.631444,2.758610,5.821964", "1.568310,6.077213,2.842349"],
    #         ["map2.txt", "1.716350,6.013881,1.933095", "1.793361,3.020136,0.962692"],
    #         ["map2.txt", "1.465272,6.219889,3.515330", "1.103914,1.781948,3.520135"],
    #         ["map2.txt", "0.639980,2.697196,2.725853", "0.333075,2.136030,2.201600"],
    #         ["map2.txt", "1.773987,6.195007,1.421362", "0.145463,3.013068,0.469323"],
    #         ["map2.txt", "0.384577,2.392987,5.896695", "1.563331,1.013648,1.295301"],
    #         ["map2.txt", "0.979561,5.592707,0.303411", "0.974934,1.121529,2.354479"],
    #         ["map2.txt", "0.608812,0.337213,1.227132", "1.100912,0.910572,0.615567"],
    #         ["map2.txt", "0.473715,2.202826,4.510053", "0.212834,1.556162,4.833221"],
    #         ["map2.txt", "0.747484,1.354368,1.344332", "0.066936,1.344987,5.744001"],
    #         ["map2.txt", "1.642982,1.129584,4.492541", "1.292501,2.329151,1.976959"],
    #         ["map2.txt", "1.377602,3.132780,3.062537", "1.674962,2.762285,6.142842"],
    #         ["map2.txt", "0.297432,3.340720,1.550566", "1.160818,3.975083,1.416427"],
    #         ]

    scores = []
    for aPlanner in [0, 1, 2, 3]:
        print("\nTESTING " + plannerList[aPlanner] + "\n")
        for i, data in enumerate(problems):
            inputMap, startPos, goalPos = [*data]
            numDOFs = len(startPos.split(","))
            outputSolutionFile = "../output/grader_out/tmp.txt"
            commandPlan = "{} {} {} {} {} {} {}".format(
                executablePath,
                inputMap, numDOFs, startPos, goalPos,
                aPlanner, outputSolutionFile)
            print("EXECUTING: " + str(commandPlan))
            commandVerify = "./../build/verifier {} {} {} {} {}".format(
                inputMap, numDOFs, startPos, goalPos,
                outputSolutionFile)
            print("EXECUTING: " + str(commandVerify))
            try:
                start = timer()
                subprocess.run(commandPlan.split(" "), check=True) # True if want to see failure errors
                timespent = timer() - start
                returncode = subprocess.run(commandVerify.split(" "), check=False).returncode
                if returncode != 0:
                    print("Returned an invalid solution")
                
                ### Calculate the cost from their solution
                with open(outputSolutionFile) as f:
                    line = f.readline().rstrip()  # filepath of the map
                    solution = []
                    for line in f:
                        solution.append(line.split(",")[:-1]) # :-1 to drop trailing comma
                    solution = np.asarray(solution).astype(float)
                    numSteps = solution.shape[0]

                    ## Cost is sum of all joint angle movements
                    difsPos = np.abs(solution[1:,]-solution[:-1,])
                    cost = np.minimum(difsPos, np.abs(2*np.pi - difsPos)).sum()

                    success = returncode == 0
                    scores.append([aPlanner, inputMap, i, numSteps, cost, timespent, success])
            
                ### Visualize their results
                commandViz = "python visualizer.py ../output/grader_out/tmp.txt --gifFilepath=../output/grader_out/grader_{}{}.gif".format(plannerList[aPlanner], i)
                commandViz += " --incPrev=1"
                subprocess.run(commandViz.split(" "), check=True) # True if want to see failure errors
            except Exception as exc:
                print("Failed: {} !!".format(exc))
                scores.append([aPlanner, inputMap, i, -1, -1, timespent, False])

    ### Save all the scores into a csv to compute overall grades
    df = pd.DataFrame(scores, columns=["planner", "mapName", "problemIndex", "numSteps", "cost", "timespent", "success"])
    df.to_csv(gradingCSV, index=False)
            

if __name__ == "__main__":
    graderMain("./../build/planner", "../output/grader_out/grader_results.csv")