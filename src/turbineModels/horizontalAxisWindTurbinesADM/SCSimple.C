//Paul Fleming
//12/07/2016
//
// "SCSimple.C" a simple implementation of a super controller meant to work in Loveseat


// OVERALL IDEA: SIMPLE CONTROLLER DOES NOTHING WITH THE INPUTS EXCEPT EXTRACT THE TIME AND DETERMINE SETPOINTS FOR YAW AND PITCH

// INPUT/OUTPUT
// inputArray, not used for now
// outputArray, for each turbine, yaw angle (compass deg) and minimum pitch (deg)
// simTime, needed to know what time it is
// numTurbines, number of turbines

#include <stdio.h>
#include <stdlib.h> // for malloc

// The main function of the super controller
void SCSimple(float * inputArray, float * outputArray, float simTime,  int numTurbines)
{

	//Define variables for handling SC_Struct
	char buffer[500];

	//Define firstTimecheck and file variables
	static int isFirstCall = 1;
	FILE *fp; //File ID
	//static int numRows;
	static int numRows = 0;

	//Declare reference vectors
	static float timeVec[1000];
	static int turbineVec[1000];
	static float yawVec[1000];
	//static float IPCVec[1000];
	//static float YMVec[1000];
	//static float TMVec[1000];
	//static float TorqueScalingVec[1000];
	static float PitchMinVec[1000];
	float temp;
	//int tempInt;

	//fprintf(stderr,"%s %d start \n",__FILE__,__LINE__);

	//Initialize row and declare memory
	if (isFirstCall == 1)
	{
		fprintf(stderr,"%s %d FIRST CALL \n",__FILE__,__LINE__);

		//Open file for writing
		if((fp=fopen("SC_INPUT.txt", "r")) == NULL) {
			//fprintf(stderr,"%s %d CANT OPEN INPUT\n",__FILE__,__LINE__);
			fprintf(stderr,"Cannot open SC Input File.\n");
			//      exit(1);
		}

		//fscanf(fp, "%d", &numRows);

		//fprintf(stderr,"%s %d NUMROWS =  %d\n",__FILE__,__LINE__,numRows);

		//Advance past the first comment line
		fgets(buffer,500,fp);


		//Now read the tables until END is found 
		//for (int i = 0; i < numRows; i++)
		for (int i = 0; 1 > 0; i++)
		{
			if(fscanf (fp, "%f", &temp) < 1)
				break;
			timeVec[i] = temp;
			fscanf (fp, "%f", &temp);
			turbineVec[i] = temp;
			fscanf (fp, "%f", &temp);
			yawVec[i] = temp;
			//fscanf (fp, "%f", &temp);
			//IPCVec[i] = temp;
			//fscanf (fp, "%f", &temp);
			//YMVec[i] = temp;
			//fscanf (fp, "%f", &temp);
			//TMVec[i] = temp;
			//fscanf (fp, "%f", &temp);
			//TorqueScalingVec[i] = temp;
			fscanf (fp, "%f", &temp);
			PitchMinVec[i] = temp;
			//fprintf(stderr,"%s %d INPUTS %f %d %f %f %f %f %f %f \n",__FILE__,__LINE__,timeVec[i],turbineVec[i],yawVec[i],IPCVec[i],YMVec[i],TMVec[i],TorqueScalingVec[i], PitchMinVec[i]);
			fprintf(stderr,"%s %d INPUTS %f %d %f %f \n",__FILE__,__LINE__,timeVec[i],turbineVec[i],yawVec[i],PitchMinVec[i]);
			numRows++;	
	}

		//Close the file
		fclose(fp);

		isFirstCall = 0; //Don't repeat
	}


	//Extract information from the current calling turbine
	//turbineID = local_SC_Struct->turbineID;
	//arrayLength = local_SC_Struct->length;
	//turbineLocalTime = local_SC_Struct->localtime;
	//data = local_SC_Struct->data;

	//Now loop through the possible conditions and apply references
	for (int i = 0; i < numRows; i++)
	{
		for (int turbineID = 0; turbineID < numTurbines; turbineID ++)
		{
		//fprintf(stderr,"%s %d tID %d turbVec %d localTime %f timeVec %f \n",__FILE__,__LINE__,turbineID,turbineVec[i],simTime,timeVec[i]);

		if (turbineID == turbineVec[i] && simTime >= timeVec[i]) // This condition applies to this turbine at this time
			{
				//fprintf(stderr,"%s %d UPDATING OUTPUT...........................................\n",__FILE__,__LINE__);

				// Set the yaw setpoint
				outputArray[turbineID * 2] = yawVec[i];
				outputArray[turbineID * 2 + 1] = PitchMinVec[i];

				//local_SC_Struct->data[38] = yawVec[i];
				//local_SC_Struct->data[39] = IPCVec[i];
				//local_SC_Struct->data[40] = YMVec[i];
				//local_SC_Struct->data[41] = TMVec[i];
				//local_SC_Struct->data[42] = TorqueScalingVec[i];
				//local_SC_Struct->data[43] = PitchMinVec[i];
			}
		}
	}

	//Write data to output
	//numTurbines = local_SC_Struct->numTurbs;
	//const char * SENSOR_LIST[SENSOR_LENGTH] = {"Wind Speed (m/s)","Power (W)","GenSpeed (rad/s)","Previous Generator Torque (nM)","IPC_ON_RAMP (-)","Blade 1 Pitch (Rad)","Blade 2 Pitch (Rad)","Blade 3 Pitch (Rad)","Blade 1 OOP Bending (Nm)","Blade 2 OOP Bending (Nm)","Blade 3 OOP Bending (Nm)","LSS Torque (Nm)","Tower Fore-Aft Accel (m/s2)","Tower Side-Side Accel (m/s2)","Rotating Hub My (Nm)","Rotating Hub Mz (Nm)","Fixed Hub My (Nm)","Fixed Hub Mz (Nm)","Yaw Bearing My (Nm)","Yaw Bearing Mz (Nm)","Yaw Angle (Degrees)","Yaw Error (Degrees)","Myd1p (Nm)","Myq1p (Nm)","Myd2p (Nm)","Myq2p (Nm)","d1Int (rad)","q1Int (rad)","d2Int (rad)","q2Int (rad)","RootMyb1P1Notch","RootMyb2P1Notch","RootMyb3P1Notch","Blade 1 In-Plane Bending (Nm)","Blade 2 In-Plane Bending (Nm)","Blade 3 In-Plane Bending (Nm)","Fore-aft Tower Bending (Nm)","Side-side Tower Bending (Nm)","SC YAW Ref (Degrees)","SC_IPC_ON","YawMomRef (Nm)","TiltMomRef (Nm)","TorqueScaling (-)","Minimum Pitch Ref (Rad)"};
	//writeTurbineData(local_SC_Struct,numTurbines,SENSOR_LIST);//Write out the turbine data

}
