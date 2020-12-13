# Pynamics_demo
Pynamics Demo Repository for RA-L submission. 
## Install Pynamics
* go to https://github.com/idealabasu/code_pynamics, download the code;
* install Pynamics via python setup.py install
* Test Pynamics in your Python IDE using the following command:
```python
    pynamics.__init__("Test")
 ```
## Install FFmpeg
    In Pynamics, we use FFmpeg to generate the simulation video. Please install FFmpeg to see the output video.
 ## Install PyCMA
  Please follow the insturction on: https://github.com/CMA-ES/pycma
## Using a Python IDE
    We commancd using Spyder for running the demo codes. Spyder can be downloaded from: https://www.spyder-ide.org/

## Parameter fitting: Use experimental data to find optimized damping ratio (Section IV-C)
    Please run: use_exp_fit_plate_dragging.py  In Spyder, after open the .py file, Run file or press F5 on keyboard to run this simulation.

## Simulation: Origami plate dragging back and forth(Large Amplitude in Section IV-B)
    Please run plate_dragging_demo_origami_large.py.

## Simulation: Soft plate dragging back and forth(Large Amplitude in Section IV-B)
    Please run: plate_dragging_demo_soft_large_amp.py

## Simulation: Soft plate dragging back and forth(Small Amplitude in Section IV-B)
    Please run: plate_dragging_demo_soft_small_amp.py

## Simulation: Robot swimming example (Section IV-D)
    Please Run: robot_swimming_demo.py

## Maxmizing Swimming efficiency (Section IV-D)
    Please Run: optimize_robot.py
