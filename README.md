# Sampling_based_planners_for_robot_arm
HW to plan path for mult dof arm from start to goal using various sampling based palnning alorithms, RRT, RRT-Connect, RRT* and PRM for course 16-782 "Planning &amp; Decision-making in Robotics"

## How to Run the Code

1. **Build the Executable**:
   - Open a terminal and navigate to the code directory.
   - Create a `build` directory:
     ```bash
     mkdir build
     cd build
     ```
   - Use CMake to build the executable:
     ```bash
     cmake ..
     cmake --build . --config Release
     ```
     Alternatively, you can compile using g++:
     ```bash
     g++ ../src/planner.cpp -o planner
     ```

2. **Run the Planner**:
   - Use the following command to execute the planner:
     ```bash
     ./planner [mapFile] [numofDOFs] [startAnglesCommaSeparated] [goalAnglesCommaSeparated] [whichPlanner] [outputFile]
     ```
   - Example:
     ```bash
     ./planner map1.txt 5 1.57,0.78,1.57,0.78,1.57 0.392,2.35,3.14,2.82,4.71 2 myOutput.txt
     ```

3. **Visualize the Output**:
   - Use the provided Python script to generate a visualization of the planned path:
     ```bash
     python scripts/visualizer.py myOutput.txt --gifFilepath myGif.gif
     ```
   - This will create a GIF file (`myGif.gif`) visualizing the motion plan.

