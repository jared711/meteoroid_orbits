# MODA
### Meteor Orbit Determination Application
MODA is a MATLAB application and accompanying library that was designed to allow meteor scientists to quickly and accurately determine the orbital elements of an observed meteoroid. It allows the user to input data observed in different coordinates and reference frames and see the effects of various perturbations including third-body effects, atmospheric drag, solar radiation pressure (SRP), and Earth oblateness effects (J2).

## Requirements
MATLAB 9+

## Installation
You can download the git repository through the terminal
```shell
git clone https://github.com/jared711/moda.git
```
Or you can download the code as a zip file and extract it to your directory of choice.

### Run in MATLAB
Navigate to the moda/ directory and run the MATLAB command
```shell
moda
```

# Usage
We have included a series of pre-programmed meteoroids that will auto-populate the required initial conditions, time data, and radar parameters. Try selecting a few of these to see how the parameters change.

## Inputs

### Frame
Choose between the East-North-Up (ENU), Earth-Centered Earth-Fixed (ECEF), or Earth-Centered Inertial (ECI) frames.

### Coordinates
Choose between cartesian, spherical, or lat/lon/az/el coordinates.

### Initial Conditions
The labels for each initial condition box will change to fit the coordinates type selected. Enter the value of each state in the left column and the standard deviation in the right column. The units of measure are displayed in the labels. If no standard deviation is known, enter 0.

### Time
Enter the time of observation of the meteoroid. The user can select the date from a calendar widget or directly enter the date in the YYYY-MM-DD format.

### Radar
The location of the radar can be input manually or autopopulated by selecting a common radar from the drop-down list.

### Mass and Area
The mass and reference area of the meteoroid are optional parameters that are only required if the Drag and SRP perturbations are turned on.

### Perturbations
Users can turn on various perturbation effects, including Three-Body, atmospheric drag, solar radiation pressure (SRP), and Earth oblateness effects (J2).

## Outputs

### Orbital Elements
The orbital elements and corresponding standard deviations are shown on the right side of the plot. The user can change the amount of significant digits used to display these values with the textbox beneath the orbital elements

### Monte Carlo
The default error propagation method uses a linearizing approximation to compute the errors at each step of the process. This approximation is good for small error values but can become less accurate for larger standard deviation values. As a backup, users can select the Monte Carlo method and change the number of particles used (N). This will randomly generate a set of N meteoroid initial conditions based on the initial state and standard deviations, compute the orbital elements for each of these initial conditions, and compute the standard deviations from the population of orbital elements. This method will take significantly longer to run than the default covariance transform method.

### Messages
Warning messages and useful information, such as the computation time, are shown in the messages box.
