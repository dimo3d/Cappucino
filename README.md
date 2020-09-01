
This is the source code for the reference implementation for the paper "Efficient 2D Simulation on Moving 3D Surfaces".

(Morgenroth et al. 2020)

The code was tested on Windows with Houdini 17.5.x and 18.0.x. Binaries of our Plugin for Houdini 17.5 and 18 on Windows 10 can be found here:
https://doi.org/10.5281/zenodo.4009208


This also works with the free apprentice version of Houdini. 
This can be downloaded at https://www.sidefx.com/tutorials/houdini-download-and-install-windows/

Licence
-------
This code is provided under MIT License.

Build procedure
---------------
Tested with Visual Studio 2017
Use the provided visual studio project files. The cmake files are not yet updated. 

Plugin-installation
-------------------
Simply copy the according dll (Cappucino.dll) into your houdini dso folder.
Typically located somewhere here:
C:\Program Files\Side Effects Software\Houdini 1x.x.xxx\houdini\dso


Usage
-----
Open one of the following files:
buoyancy_flow.hipnc
rotating_sphere.hipnc

You will find one Geometry node inside these files.
On this geometry node, we added a parameter tab called SimAttributes.
In this parameter tab the most important simulation attributes can be manipulated.

For both scenarios we set the parameter as presented in the paper, except the simulation resolution in the rotating_sphere example. 
There we chose a coarser one, for faster testing. You can adjust the resulution with the "simresolution" attribute.

Parameter List
--------------

Fluid attributes:
Vorticity			Changes the fluids vorticity. Note that this parameter is resolution dependent. 
				When changing the sim resolution you need to adjust this parameter aswell.
Velocity Diffusion		The fluids dynamic viscosity constant.
Coupling Friction		Coupling parameter s_1 from the paper. Enlarges adhesion forces
Coupling Curl			Coupling parameter s_2 from the paper. Controls the strength of the Jacobian based coupling
Coupling Noise			With this parameter you can add some noise to the coupling, introducing more noise in the velocity field and
				therefore more fine scale details on the simulation.

Buoyancy attributes:
Gravity 			Enables gravity to be scaled as described in the Boussinesq approximation.
Heating Rate 			Strength of the temperature transfer from the outer constraints for the hot constraint.
Cooling Rate 			Strength of the temperature transfer from the outer constraints for the cold constraint.
Heat Dissipation		Temperature change due to heat induced velocity divergence.
Heat Diffusion			Reduces relative temperature differences.
Cool Temperature  		Lowest allowed temperature and temperature of the outer cooling constraint
High Temperature 		Highest allowed temperature and temperature of the outer heat constraint
Reference Temperature   	t_0 of the Boussinesq approximation. In areas of this temperature basically 0 gravity is present.

General Attributes:
simresolution			The VDB-Grid Resolution.
Interpolation Method		The CPM-Interpolation Method
CPM Cells			Width of the Narrowband
Advection Scheme		Integrator used for the Advection of the VDB fields
Advection Substeps		Number of advection substeps (Should in general be 1, for more details rather increase the substeps attribute)
Substeps 			Number of Timesteps per Frame (Together with Houdini's $FPS attribute, it determines the simulation step size)
Reset simulation 		Resets the simulation
