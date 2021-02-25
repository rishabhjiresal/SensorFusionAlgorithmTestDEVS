if [ -z $1 ]
then
  selection="Please run the script with args. After ./RUNME.sh, type 'sim' for simulation, 'emb' for project compilation for real time hardware or type 'embf' to compile and flash for realtime hardware."
elif [ -n $1 ]
then
  selection=$1
fi
case $selection in
   "sim") echo "Running the simulation : ";cd top_model; make clean; make all; ./Testing_Algorithm_TOP;echo "Please check the top_model folder for output";;
   "emb") echo "Building for real time hardware: ";mbed new .;cd top_model; make clean; make embedded;;
   "embf") echo "Building for real time hardware Flashing will be done after building, please connect your microcontroller board ";mbed new .;cd top_model; make clean; make embedded; make flash;;
   *) echo "$selection";;
esac
