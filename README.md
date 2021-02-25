### Sensor Fusion Algorithm Test DEVS ###


### REAL TIME CADMIUM INSTALL ###

Clone this repo into an empty folder, ensure you download the MBed submodules:

> git submodule update --init --recursive

Run this to install dependencies

> SensorFusionAlgorithmTestDEVS/install.sh
### RUNNING THE MODEL USING AUTOMATED SCRIPT ###

First, open the terminal in the SensorFusionAlgorithmTestDEVS folder. Type
> chmod +x install.sh

and then type

> sudo ./install.sh

This is very important since all the dependencies and requirements should be installed properly to be able to compile the project successfully. (Note: Please extract the archive to the home directory for smooth experience)

There is a script written for compiling, simulating and flashing the project for the NUCLEO Board. Open terminal in the folder SensorFusionAlgorithmTestDEVS and type

> chmod +x RUNME.sh.

To build and run the simulation, type

> ./RUNME.sh sim

This will build and run the simulation and the output can be seen in the top_model folder in outputs and log will be generated in top_model folder.

To build the project for NUCLEO board, type

> ./RUNME.sh emb (Note: This might take some time)

To build and flash the project on the NUCLEO Board, type

> ./RUNME.sh embf

(Note: If you have a different board other than the F401RE or F446RE please make sure you change the COMPILE_TARGET and FLASH_TARGET in the makefile in top_model folder.)

### SIMULATE MODELS ###

> cd SensorFusionAlgorithmTestDEVS/top_model/

> make clean; make all

> ./Testing_Algorithm_TOP

This will run the standard Cadmium simulator. Cadmium logs will be generated in SensorFusionAlgorithmTestDEVS/top_model/SensorFusion_Cadmium_output.txt
The pin's inputs are stored in SensorFusionAlgorithmTestDEVS/top_model/inputs. The value of the output pins will be in SensorFusionAlgorithmTestDEVS/top_model/inputs.
SVEC (Simulation Visualization for Embedded Cadmium) is a python GUI that parses these files and steps through the simulation to help debug the models.


### RUN MODELS ON TARGET PLATFORM ###

If your target platform *is not* the Nucleo-STM32F401, you will need to change the COMPILE_TARGET / FLASH_TARGET in the make file.

> cd SensorFusionAlgorithmTestDEVS

> mbed new .

> cd SensorFusionAlgorithmTestDEVS/top_model/

> make clean; make embedded; make flash;
