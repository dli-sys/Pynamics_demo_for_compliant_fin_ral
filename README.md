# Pynamics_demo
Pynamics Demo Repository for RA-L submission. 
## Install Pynamics
* Go to https://github.com/idealabasu/code_pynamics, download the code;
* Install Pynamics via python setup.py install. We recommanded this method to get latest version of Pynamics;
* In Pynamics, we use FFmpeg to generate the simulaticon video. Please install FFmpeg to see the output video. Official website: https://ffmpeg.org/
* To see if Pynamics is correctly installed, type 
```python
import pynamics
```
If no error or exception was given, Pynamics is correctly installed.

 ## Install PyCMA
Please follow the insturction on: https://github.com/CMA-ES/pycma
## Using a Python IDE
We commancd using Spyder for running the demo code. Our code is only tested in this environment. Spyder can be downloaded from: https://www.spyder-ide.org/
## Before run the demos:
Download this repository and unzip, please make sure fit_qs.py is in the same directory, which will be required in later simulation
## Parameter fitting: Use experimental data to find optimized damping ratio (Section IV-C)
Please run: use_exp_fit_plate_dragging.py  In Spyder, after open the .py file, Run file or press F5 on keyboard to run this simulation.
fit_qs.py is the file where we store our experimental data, please make sure it's placed in the same directory.

## Simulation: Origami plate dragging back and forth(Large Amplitude in Section IV-B)
Please run plate_dragging_demo_origami_large.py
Demo video:
![alt text](https://github.com/gdbbzq/Pynamics_demo/blob/main/demo_videos/plate_dragging_demo_origami_large.gif)

## Simulation: Soft plate dragging back and forth(Large Amplitude in Section IV-B)
Please run: plate_dragging_demo_soft_large_amp.py
Demo video:
![alt text](https://github.com/gdbbzq/Pynamics_demo/blob/main/demo_videos/plate_dragging_demo_soft_large_amp.gif)

## Simulation: Soft plate dragging back and forth(Small Amplitude in Section IV-B)
Please run: plate_dragging_demo_soft_small_amp.py
Webpage preview might not show the video correctly, we recommand download the video. Demo video:
![alt text](https://github.com/gdbbzq/Pynamics_demo/blob/main/demo_videos/plate_dragging_demo_soft_small_amp.gif)

## Simulation: Robot swimming example (Section IV-D)
Please Run: robot_swimming_demo.py
Webpage preview might not show the video correctly, we recommand download the video. Demo video:
![alt text](https://github.com/gdbbzq/Pynamics_demo/blob/main/demo_videos/robot_swimming_demo.gif)

## Maxmizing Swimming efficiency (Section IV-D)
Please Run: optimize_robot.py
