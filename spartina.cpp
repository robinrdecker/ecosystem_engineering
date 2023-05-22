/* Robin Decker
12/14/2017

Spartina

Numerically analyze Ecosystem Engineering model

The header files "nr3.h", "interp_1d.h", romberg.h", "quadrature.h" are propriety software from Numerical Recipes, 3rd edition.
If you don't have access to this software, you can replace any qromb statements with your personal function that performs
romberg integration, or any other quadrature technique that you prefer.
*/

#include "Rcpp.h" //for  R

#include "nr3.h" //numerical recipes data structures
#include "interp_1d.h" //interpolation methods
#include "romberg.h" //Romberg integration
#include "quadrature.h" //base numerical integration
#include <math.h> //power function, exponential function
#include <iostream> //print to console

#include <unistd.h>

using namespace Rcpp;

//Global parameters:

//The initial population level
const double INITIALPOP = 0.01;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Function declarations:

//Function that does all of the work of the program
List runModel(const int &subintervals, const int &timesteps, const double &d, const double &m, const double &l,
	const double &mu, const double &sigma, const double &rMax, const double &b, const double &c, 
	const int &initialHabitatQuality, const bool &doEngineer, const bool &doErosion, const bool &doriseSeaLevel,
	const bool &doGrow);

//Sets up initial u
void setUInitial(VecDoub &uInitial, const int &subintervals, const int &initialHabitatQuality, 
	const VecDoub &habitatQuality);

//Sets up initial o
void setOInitial(VecDoub &oInitial, const int &subintervals, const VecDoub &habitatQuality);

//Performs the erosion step
void erode(const VecDoub &uBefore, const VecDoub &oBefore, VecDoub &uAfter, VecDoub &oAfter, const double &d,
	const int &subintervals, const VecDoub &habitatQuality, const VecDoub &uBeforeAtoms, const VecDoub &oBeforeAtoms,
	VecDoub &uAfterAtoms, VecDoub &oAfterAtoms, const VecDoub &habitatQualityAtoms);

//Performs the engineering step
void engineer(const VecDoub &uBefore, const VecDoub &oBefore, VecDoub &uAfter, VecDoub &oAfter,
	const double &m, const double &l, const int &subintervals, const VecDoub &habitatQuality,
	const VecDoub &uBeforeAtoms, const VecDoub &oBeforeAtoms,
	VecDoub &uAfterAtoms, VecDoub &oAfterAtoms, const VecDoub &habitatQualityAtoms);

//Performs the sea-level rise step
void riseSeaLevel(const VecDoub &uBefore, const VecDoub &oBefore, VecDoub &uAfter, VecDoub &oAfter,
	const double &c, const int &subintervals, const VecDoub &habitatQuality,  const VecDoub &uBeforeAtoms, const VecDoub &oBeforeAtoms,
	VecDoub &uAfterAtoms, VecDoub &oAfterAtoms, VecDoub &habitatQualityAtoms, const int &timestep);

//Performs the growth step
void grow(const VecDoub &uBefore, const VecDoub &oBefore, VecDoub &uAfter, VecDoub &oAfter, const double &mu,
	const double &sigma, const double &rMax, const double &b, const int &subintervals, const VecDoub &habitatQuality,
	const VecDoub &uBeforeAtoms, const VecDoub &oBeforeAtoms,
	VecDoub &uAfterAtoms, VecDoub &oAfterAtoms, const VecDoub &habitatQualityAtoms);

void disperse(const VecDoub &oBeforeAtoms, const VecDoub &uBeforeAtoms, const VecDoub &habitatQuality, const VecDoub &habitatQualityAtoms, 
	VecDoub &AtomsAfter, const int &indexOfSource, const double &b, const double &mu, const double &sigma, const double &rMax,
	const VecDoub &oBefore, const VecDoub &uBefore, VecDoub &oAfter, const VecDoub &densityDependenceVector);

double getDensityDependence(const VecDoub &oBefore, const VecDoub &uBefore, const VecDoub &oBeforeAtoms, const VecDoub &uBeforeAtoms,
    const double &start, double &end, const VecDoub &habitatQuality, const VecDoub &habitatQualityAtoms);

void setDensityDependenceVector(const VecDoub &oBefore, const VecDoub &uBefore, const VecDoub &oBeforeAtoms, const VecDoub &uBeforeAtoms, 
	const VecDoub &habitatQuality, const VecDoub &habitatQualityAtoms, VecDoub &densityDependenceVector);

//A function describing the dependence of habitat modification on population density, used in engineering step
double h(const VecDoub &oBefore, const double &l, const double &sum, const int &subintervals, const double &hq, //habitat quality value
	const VecDoub &oBeforeAtoms, const VecDoub &habitatQualityAtoms, const VecDoub &habitatQuality);

//Function describing the growth rate, which is needed for the population growth function, used in the growth step
double r(const double &x, const double &mu, const double &sigma, const double &rMax);

//Function that performs interpolation
double interpolate(const VecDoub &X, const VecDoub &Y, const double x);

//Converts NumericVector types to VecDoub types
VecDoub NumericVectorToVecDoub(const NumericVector &vec);

//Converts VecDoub types to NumericVector types
NumericVector VecDoubToNumericVector(const VecDoub &vec);

//Function that records each iteration to the NumericVector
void saveToRVector(const VecDoub &cStructure, NumericVector &rStructure, const int &pointer);

//Function to fill out t
void setT(const int &subintervals, const int &timesteps, NumericVector &t);

//Function that fills in habitat quality vector
void setHabitatQuality(VecDoub &habitatQuality, const int &subintervals);

void setX(const int &subintervals, const int &timesteps, NumericVector &x, const VecDoub &habitatQuality);

//Functions for atom data structures
void setHabitatQualityAtoms(VecDoub &habitatQualityAtoms);
void setUInitialAtoms(VecDoub &uInitialAtoms);
void setOInitialAtoms(VecDoub &oInitialAtoms);

void setTAtoms(NumericVector &tAtoms, const int &timesteps, const double &c);
void setXAtoms(NumericVector &xAtoms, const int &timesteps, const double &c);

void erodeLeft(const VecDoub &AtomsInitial, const VecDoub &habitatQuality, const VecDoub &habitatQualityAtoms, VecDoub &Eroded,
	VecDoub &AtomsEroded, const int &indexOfSource, const double &d);

void engineerRight(const VecDoub &BeforeAtoms, const VecDoub &habitatQuality, const VecDoub &habitatQualityAtoms, VecDoub &After,
	VecDoub &AtomsAfter, const int &indexOfSource, const double &b);

//Finds the number of atoms based on the number of timesteps run and the speed of sea level rise
int getNumOfAtoms(const int &timesteps, const double &c);

double getTotalOccupiedArea(const VecDoub &oBefore, const VecDoub &oBeforeAtoms, const VecDoub &habitatQuality);

double signum(double &x);

void finalize(VecDoub &before);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Functor definitions:

//Simple integrand
struct simpleIntegrand {
	simpleIntegrand(const VecDoub &habitatQuality, const VecDoub &dependentVar) :
		habitatQuality(habitatQuality), dependentVar(dependentVar) {}
	double operator()(double x) const {
		return interpolate(habitatQuality, dependentVar, x);
}

private:
	VecDoub habitatQuality;
	VecDoub dependentVar;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*Function that does all of the work of the program*/

// [[Rcpp::export()]]
List runModel(const int &subintervals, const int &timesteps, const double &d, const double &m, const double &l,
	const double &mu, const double &sigma, const double &rMax, const double &b, const double &c, 
	const int &initialHabitatQuality, const bool &doEngineer, const bool &doErosion, const bool &doRiseSeaLevel,
	const bool &doGrow){
	
	//Find the number of atoms
	int numOfAtoms = getNumOfAtoms(timesteps, c); //number of atoms summed across all timesteps

	//Create C++ vectors

	//Create a vector to store values of habitat quality or "sediment level" and initialize this vector
	VecDoub habitatQuality(subintervals);
	setHabitatQuality(habitatQuality, subintervals);

	//Create a vector to store values of habitat quality for the finite number of atoms, and initialize this vector
	VecDoub habitatQualityAtoms(2); //2 is the number of initial atoms (one at -2 and one at 1)
	setHabitatQualityAtoms(habitatQualityAtoms);
	
	//Create vectors to store the initial distribution of unoccupied area and occupied area at each tidal height, and
	//initialize both of these vectors
	VecDoub uBefore(subintervals);
	setUInitial(uBefore, subintervals, initialHabitatQuality, habitatQuality);
	VecDoub oBefore(subintervals);
	setOInitial(oBefore, subintervals, habitatQuality);

	//Create vectors to store the initial vales of unoccupied and occupied area at each atom, and
	//initialize both of these vectors
	VecDoub uBeforeAtoms(2); //2 is the number of initial atoms (one at -2 and one at 1)
	setUInitialAtoms(uBeforeAtoms);
	VecDoub oBeforeAtoms(2); //2 is the number of initial atoms (one at -2 and one at 1)
	setOInitialAtoms(oBeforeAtoms);

	//After the growth step, we store the new area unoccupied and occupied at each tidal height in these vectors
	//There are 2 vectors for the distribution of values, and 2 vectors for the atoms
	VecDoub uAfter(subintervals);
	uAfter = uBefore; //In case none of the 4 steps are done, these "after" data structures still need to be set
	VecDoub oAfter(subintervals);
	oAfter = oBefore;
	VecDoub uAfterAtoms(2); //2 is the number of initial atoms (one at -2 and one at 1)
	uAfterAtoms = uBeforeAtoms;
	VecDoub oAfterAtoms(2); //2 is the number of initial atoms (one at -2 and one at 1)
	oAfterAtoms = oBeforeAtoms;

	//Declare R vectors
	NumericVector u(subintervals*(timesteps + 1));
	NumericVector o(subintervals*(timesteps + 1));
	NumericVector times(subintervals*(timesteps + 1));
	NumericVector x(subintervals*(timesteps + 1)); 

	//Declare R vectors for the atoms
	NumericVector uAtoms(numOfAtoms); //numOfAtoms is the total number of atoms summed across all timesteps
	NumericVector oAtoms(numOfAtoms);
	NumericVector tAtoms(numOfAtoms);
	NumericVector xAtoms(numOfAtoms);

	//Set t
	setT(subintervals, timesteps, times);
	setX(subintervals, timesteps, x, habitatQuality);


	setTAtoms(tAtoms, timesteps, c);
	setXAtoms(xAtoms, timesteps, c);

	int pointer = 0;

	//Find u and o, one timestep at a time:
	for (int i = 0; i < timesteps; i++) {

		cout << "Timestep " << i << endl; //Keep track of progress

		//Save the initial data structures (not the atoms)
		saveToRVector(uBefore, u, i*subintervals);
		saveToRVector(oBefore, o, i*subintervals);

		//Save the initial atoms
		saveToRVector(uBeforeAtoms, uAtoms, pointer);
		saveToRVector(oBeforeAtoms, oAtoms, pointer);
		
		pointer += uBeforeAtoms.size(); 

		//saveToRVector(habitatQuality, x, subintervals, i); Need to do this at the same time as the t vector for R

		//Perform all 4 steps
		if (doErosion){

			erode(uBefore, oBefore, uAfter, oAfter, d, subintervals, habitatQuality, uBeforeAtoms, oBeforeAtoms,
				uAfterAtoms, oAfterAtoms, habitatQualityAtoms);
			finalize(uAfter);
			finalize(oAfter);
			finalize(uAfterAtoms);
			finalize(oAfterAtoms);
			cout << "Erosion finished" << endl;

			//The vectors at the last step become the new vectors for the next timestep
			uBefore = uAfter;
			oBefore = oAfter;
			uBeforeAtoms = uAfterAtoms;
			oBeforeAtoms = oAfterAtoms;
		}
		if (doEngineer){

			engineer(uBefore, oBefore, uAfter, oAfter, m, l, subintervals, habitatQuality,
				uBeforeAtoms, oBeforeAtoms, uAfterAtoms, oAfterAtoms, habitatQualityAtoms);
			finalize(uAfter);
			finalize(oAfter);
			finalize(uAfterAtoms);
			finalize(oAfterAtoms);
			cout << "Engineering finished" << endl;

			//The vectors at the last step become the new vectors for the next timestep
			uBefore = uAfter;
			oBefore = oAfter;
			uBeforeAtoms = uAfterAtoms;
			oBeforeAtoms = oAfterAtoms;
		}
		if (doRiseSeaLevel){

			riseSeaLevel(uBefore, oBefore, uAfter, oAfter, c, subintervals, habitatQuality,  uBeforeAtoms, oBeforeAtoms,
				uAfterAtoms, oAfterAtoms, habitatQualityAtoms, i);
			finalize(uAfter);
			finalize(oAfter);
			finalize(uAfterAtoms);
			finalize(oAfterAtoms);
			
			cout << "Sea level rise finished" << endl;

			//The vectors at the last step become the new vectors for the next timestep
			uBefore = uAfter;
			oBefore = oAfter;
			uBeforeAtoms = uAfterAtoms;
			oBeforeAtoms = oAfterAtoms;
		}
		if (doGrow){

			grow(uBefore, oBefore, uAfter, oAfter, mu, sigma, rMax, b, subintervals, habitatQuality,
				uBeforeAtoms, oBeforeAtoms, uAfterAtoms, oAfterAtoms, habitatQualityAtoms);
			finalize(uAfter);
			finalize(oAfter);
			finalize(uAfterAtoms);
			finalize(oAfterAtoms);

			cout << "Growth finished" << endl;

			//The vectors at the last step become the new vectors for the next timestep
			uBefore = uAfter;
			oBefore = oAfter;
			uBeforeAtoms = uAfterAtoms;
			oBeforeAtoms = oAfterAtoms;
		}
	}

	//Save the last iteration (first the non-atoms, then the atoms)
	saveToRVector(uAfter, u, subintervals*timesteps);
	saveToRVector(oAfter, o, subintervals*timesteps);

	saveToRVector(uAfterAtoms, uAtoms, pointer);
	saveToRVector(oAfterAtoms, oAtoms, pointer);

	//Return the data structure to R
	DataFrame atoms = DataFrame::create(_["Timestep"] = tAtoms, _["Quality"] = xAtoms, _["Unoccupied"] = uAtoms, _["Occupied"] = oAtoms);
	DataFrame df = DataFrame::create(_["Timestep"] = times, _["Quality"] = x, _["Unoccupied"] = u, _["Occupied"] = o);
	cout << "Returning both dataframes" << endl;
	return Rcpp::List::create(_["Atoms"] = atoms, _["df"] = df);
}

//Sets up initial u
void setUInitial(VecDoub &uInitial, const int &subintervals, const int &initialHabitatQuality, 
	const VecDoub &habitatQuality) {

	//Initially neutral
	if (initialHabitatQuality == 1) {
		for (int i = 0; i < subintervals; i++) {
			uInitial[i] = 1.0;
		}
	} 
	else { //Initially suboptimal
		double mu;
		if (initialHabitatQuality == 0) {
			mu = -1;

		} //Initially optimal
		else { //if initialHabitatQuality == 2 or anything else
			mu = 0;
		}
		//These values make everything integrate to 4 (same area as neutral case)
		double sigma = 0.15;
		double scale = 1.8;
		double intercept = 0.55;
		for (int i = 0; i < subintervals; i++) {
			uInitial[i] = intercept + scale*(1 / pow(2 * PI*pow(sigma, 2.0), 0.5))*exp(-pow(habitatQuality[i] - mu, 2) / (2 * pow(sigma, 2.0)));
		}
	}

}

//Sets up initial o
void setOInitial(VecDoub &oInitial, const int &subintervals, const VecDoub &habitatQuality) {

	for (int i = 0; i < subintervals; i++) {
		if (habitatQuality[i] >= -1.05 & habitatQuality[i] <= -0.95){
			oInitial[i] = INITIALPOP;
		} else {
			oInitial[i] = 0.0;
		}
	}

	//oInitial[(3*subintervals)/8] = INITIALPOP;
}

//To move area out of atoms and into other atoms or the background surrounding atoms 
//(this is for one atom at a time, indexed at indexOfSource)
void erodeLeft(const VecDoub &AtomsInitial, const VecDoub &habitatQuality, const VecDoub &habitatQualityAtoms, VecDoub &Eroded, 
	VecDoub &AtomsEroded, const int &indexOfSource, const double &d) {

	double modifier = .50;

	double progeny = AtomsInitial[indexOfSource];
	double locationOfSource = habitatQualityAtoms[indexOfSource];
	int subintervals = habitatQuality.size();

	//first take care of the background
	for (int dest = subintervals - 1; dest >= 0; dest--) { //go backwards to keep track of endpoints
		if (habitatQuality[dest] < locationOfSource) { //background points must be to the left of the source atom 
			
            Eroded[dest] += modifier*progeny*(1 / d)*exp(-abs(habitatQuality[dest] - locationOfSource) / d);
         
		} //end the if statement that checks that background points are to the left of the source atom
	} //end the for loop that goes backwards to keep track of endpoints

    AtomsEroded[indexOfSource] -= modifier*progeny;
	//If we are accounting for all of the 10 percent of area lost from the atom, then we may as well just take 10 percent away from
	//the atom in one step, instead of subtracting little portions along the way
}


void engineerRight(const VecDoub &BeforeAtoms, const VecDoub &habitatQuality, const VecDoub &habitatQualityAtoms, VecDoub &After,
	VecDoub &AtomsAfter, const int &indexOfSource, const double &b){

		double modifier = .50;

		double progeny = BeforeAtoms[indexOfSource]; 
		double locationOfSource = habitatQualityAtoms[indexOfSource];
		int subintervals = habitatQuality.size();

		//first take care of the background
		if (b != 0){

		for (int dest = 0; dest < subintervals; dest++){ //go forwards to keep track of endpoints
			if (habitatQuality[dest] > locationOfSource && habitatQuality[dest] <= 1){ //background points must be greater than source
			//Also need to make sure that the destination is <= 1, since no engineering past 1

                After[dest] += modifier*progeny*(1/b)*exp(-(habitatQuality[dest]-locationOfSource)/b);

			} // end the if statement that check that the backgrund points are right of the source and less than 1

		} // end the for loop that goes forwards to keep track of endpoints

		} else {
			AtomsAfter[indexOfSource] += modifier*progeny;
		}

        AtomsAfter[indexOfSource] -= modifier*progeny;

	}


//Performs the erosion step, size of uAfterAtoms, oAfterAtoms and habitatQualityAtoms do no change
void erode(const VecDoub &uBefore, const VecDoub &oBefore, VecDoub &uAfter, VecDoub &oAfter, const double &d,
	const int &subintervals, const VecDoub &habitatQuality, const VecDoub &uBeforeAtoms, const VecDoub &oBeforeAtoms,
	VecDoub &uAfterAtoms, VecDoub &oAfterAtoms, const VecDoub &habitatQualityAtoms) {

	if (d > 0) { //If there is some erosion

		//1.) Reset uAfter and oAfter so that they are empty:
		for (int i = 0; i < subintervals; i++) {

			uAfter[i] = 0.0;
			oAfter[i] = 0.0;
		}
		//Set AtomsEroded
		uAfterAtoms = uBeforeAtoms;
		oAfterAtoms = oBeforeAtoms;

		//2.) Create the integrand (this does not include the kernel)
		double uprogeny, oprogeny, lower, upper, transferAmount, upperdiff, lowerdiff;

		//3.) Go one interval at a time, calculating the number of progeny in the interval
		for (int source = 0; source < subintervals; source++) { //If source is 0, do nothing, since everything stays there 

			uprogeny = uBefore[source];
			oprogeny = oBefore[source];

			//4.) Calculate how many of those progeny go to each of the k other intervals (dispersal kernel)
			for (int dest = 0; dest < subintervals; dest++) {
				if (habitatQuality[dest] >= -2.0 && habitatQuality[dest] < habitatQuality[source]) { //All the normal cases

                    //For centered subintervals
			        lower = habitatQuality[dest] - 0.5*(4.0 / subintervals);
			        upper = habitatQuality[dest] + 0.5*(4.0 / subintervals);

                    upperdiff = upper - habitatQuality[source];
			        lowerdiff = lower - habitatQuality[source];

                    transferAmount = exp(upperdiff/d) - exp(lowerdiff/d);

                    uAfter[dest] += uprogeny*transferAmount;
                    oAfter[dest] += oprogeny*transferAmount;
				}
			}
		} //End of for loop going source by source
        
            //4.5) For each source, some goes to the atom at -2 (at the end, we will put all leftovers in the atom)


		//5.) Finally, take care of all other atoms as the source, one at a time, using the erodeLeft method

		for (int asource = 0; asource < habitatQualityAtoms.size(); asource++) { //go atom by atom as source
			erodeLeft(oBeforeAtoms, habitatQuality, habitatQualityAtoms, oAfter, oAfterAtoms, asource, d);
			erodeLeft(uBeforeAtoms, habitatQuality, habitatQualityAtoms, uAfter, uAfterAtoms, asource, d);
		}

		//6.) Put leftovers in the atom at -2.
		uAfterAtoms[0] = 0;
		oAfterAtoms[0] = 0;

        //First, finalize
        finalize(oAfter);
        finalize(uAfter);
		finalize(uAfterAtoms);
		finalize(oAfterAtoms);

		uAfterAtoms[0] = getTotalOccupiedArea(uBefore, uBeforeAtoms, habitatQuality) - getTotalOccupiedArea(uAfter, uAfterAtoms, habitatQuality);
		oAfterAtoms[0] = getTotalOccupiedArea(oBefore, oBeforeAtoms, habitatQuality) - getTotalOccupiedArea(oAfter, oAfterAtoms, habitatQuality);

	}
	else { //if there is no erosion
		uAfter = uBefore;
		oAfter = oBefore;
		uAfterAtoms = uBeforeAtoms;
		oAfterAtoms = oBeforeAtoms;
	}
}



//Performs the engineering step, size of uAfterAtoms, oAfterAtoms and habitatQualityAtoms does not change
void engineer(const VecDoub &uBefore, const VecDoub &oBefore, VecDoub &uAfter, VecDoub &oAfter,
	const double &m, const double &l, const int &subintervals, const VecDoub &habitatQuality,
	const VecDoub &uBeforeAtoms, const VecDoub &oBeforeAtoms,
	VecDoub &uAfterAtoms, VecDoub &oAfterAtoms, const VecDoub &habitatQualityAtoms) {

		if (m > 0) { // if there is some engineering

			//1.) Reset uAfter and oAfter so that they are empty

			for (int i = 0; i < subintervals; i++){

				if (habitatQuality[i] <= 1){
					uAfter[i] = 0.0;
					oAfter[i] = 0.0;
				} else {
					uAfter[i] = uBefore[i];
					oAfter[i] = oBefore[i];
				}
			}

			//Set up initial values of AfterAtoms
			uAfterAtoms = uBeforeAtoms;
			oAfterAtoms = oBeforeAtoms;

			//2.) Create the integrand (this does not include the kernel)

			double uprogeny, oprogeny, lower, upper, transferAmount, upperdiff, lowerdiff;

			double totalArea = getTotalOccupiedArea(oBefore, oBeforeAtoms, habitatQuality); //this is the same for whole timestep
			double b; //will change for each source, but is the same for all destinations within a source

			//3.) Go one source interval at a time, calculating the number of progeny in the interval

			for (int source = 0; habitatQuality[source] <= 1; source++){ //only look at sources less than 1, as we can engineer past 1

				uprogeny = uBefore[source];
				oprogeny = oBefore[source];

				//b is the same for the entire source
				b = m*h(oBefore, l, totalArea, subintervals, habitatQuality[source], oBeforeAtoms, habitatQualityAtoms, habitatQuality);

				if (b != 0){

					//4.) Calculate what proportion of those progeny go to each of the other intevals using the dispersal kernel
					for (int dest = 0; habitatQuality[dest] <= 1; dest++){ //only look at destinations less than 1, no engineering past 1

						//First look at all of the normal cases that don't involve atoms
						if (habitatQuality[dest] > habitatQuality[source]){

                        	//For centered subintervals
				        	lower = habitatQuality[dest] - 0.5*(4.0/subintervals);
				        	upper = habitatQuality[dest] + 0.5*(4.0/subintervals);

                        	upperdiff = habitatQuality[source] - upper;
                        	lowerdiff = habitatQuality[source] - lower;

                        	transferAmount = exp(lowerdiff/b) - exp(upperdiff/b); 

							uAfter[dest] += uprogeny*transferAmount;
							oAfter[dest] += oprogeny*transferAmount;

						}
					}
				} else {
					uAfter[source] += uprogeny;
					oAfter[source] += oprogeny;
				}
				//4.5) For each source, some goes to the atom at 1
				//the atom at 1 must always be the atom with the highest habitat Quality value, because sea level rise only shifts atoms down
				//Will take care of this as leftovers in step 6

			} //end of for loop going source by source
			
			//5.) Finally, take care of all other atoms as the source, one at a time, using the engineerRight method
			for (int asource = 0; asource < habitatQualityAtoms.size(); asource++){
				double b = m*h(oBefore, l, totalArea, subintervals, habitatQualityAtoms[asource], oBeforeAtoms, habitatQualityAtoms, habitatQuality);
					engineerRight(oBeforeAtoms, habitatQuality, habitatQualityAtoms, oAfter, oAfterAtoms, asource, b);
					engineerRight(uBeforeAtoms, habitatQuality, habitatQualityAtoms, uAfter, uAfterAtoms, asource, b);

			}

			//6.) Put all of the leftovers in the atom at 1:
			int indexOfAtomOne = habitatQualityAtoms.size() - 1;
			uAfterAtoms[indexOfAtomOne] = 0;
			oAfterAtoms[indexOfAtomOne] = 0;

            //First, finalize
            finalize(oAfter);
            finalize(uAfter);
			finalize(uAfterAtoms);
			finalize(oAfterAtoms);

			double totalUBefore = getTotalOccupiedArea(uBefore, uBeforeAtoms, habitatQuality);
            double totalUAfter = getTotalOccupiedArea(uAfter, uAfterAtoms, habitatQuality);
            if (totalUBefore - totalUAfter < 0){
                uAfterAtoms[0] -= totalUAfter - totalUBefore;
                uAfterAtoms[indexOfAtomOne] = 0;
            } else {
			    uAfterAtoms[indexOfAtomOne] = totalUBefore - totalUAfter;
            }

			double totalOBefore = getTotalOccupiedArea(oBefore, oBeforeAtoms, habitatQuality);
            double totalOAfter = getTotalOccupiedArea(oAfter, oAfterAtoms, habitatQuality);
            if (totalOBefore - totalOAfter < 0){
                oAfterAtoms[0] -= totalOAfter - totalOBefore;
                oAfterAtoms[indexOfAtomOne] = 0;
            }
            else {
			    oAfterAtoms[indexOfAtomOne] = totalOBefore - totalOAfter;
            }

		} // end of "if engineering is not set to 0"
		else { //if engineering is set to 0
			uAfter = uBefore;
			oAfter = oBefore;
			uAfterAtoms = uBeforeAtoms;
			oAfterAtoms = oBeforeAtoms;
		}

	} // end function


//Performs the sea-level rise step
//This is the step where the vectors containing the atoms grow
//uAfterAtoms, oAfterAtoms and habitatQualityAtoms all grow by 1
void riseSeaLevel(const VecDoub &uBefore, const VecDoub &oBefore, VecDoub &uAfter, VecDoub &oAfter,
	const double &c, const int &subintervals, const VecDoub &habitatQuality,  const VecDoub &uBeforeAtoms, const VecDoub &oBeforeAtoms,
	VecDoub &uAfterAtoms, VecDoub &oAfterAtoms, VecDoub &habitatQualityAtoms, const int &timestep) {


	if (c > 0) { //if there is some sea level rise

		//1.) habitatQualityAtoms increases by 1 and shifts
		VecDoub newHabitatQualityAtoms(habitatQualityAtoms.size()+1); //one longer than the old vector
		newHabitatQualityAtoms[0] = -2.0;

		for(int i = newHabitatQualityAtoms.size() - 1; i > 0; i-- ){
			if (i > 1){
				newHabitatQualityAtoms[i] = habitatQualityAtoms[i-1];
			} else {
				newHabitatQualityAtoms[i] = habitatQualityAtoms[i] - c;
			}
		}

		habitatQualityAtoms = newHabitatQualityAtoms;

		//2.) newUAtoms is created, one unit longer than uBeforeAtoms, and contents are shifted

		VecDoub newUAtoms(uBeforeAtoms.size()+1);
		for (int i = 0; i < uBeforeAtoms.size(); i++){
			newUAtoms[i] = uBeforeAtoms[i];
		}
		newUAtoms[newUAtoms.size()-1] = 0;

		//3.) oAfterAtoms is created, one unit longer than oBeforeAtoms, and contents are shifted

		VecDoub newOAtoms(oBeforeAtoms.size()+1);
		for(int i = 0; i < oBeforeAtoms.size(); i++){
			newOAtoms[i] = oBeforeAtoms[i];
		}
		newOAtoms[newOAtoms.size()-1] = 0;

		//4.) The background is shifted

		for (int i = 0; i < subintervals; i++){
			if (habitatQuality[i] > 2.0-c){
				uAfter[i] = 0;
				oAfter[i] = 0;
			} else {
				uAfter[i] = interpolate(habitatQuality, uBefore, habitatQuality[i]+c);
				oAfter[i] = interpolate(habitatQuality, oBefore, habitatQuality[i]+c);
			}
		}

		//5.) The portion of the background that is shifted out of bounds is distributed to the atom at -2
		
		simpleIntegrand oIntegrand(habitatQuality, oBefore);
		simpleIntegrand uIntegrand(habitatQuality, uBefore); 

		newUAtoms[0] += qromb(uIntegrand, -2, -2 + c);
		newOAtoms[0] += qromb(oIntegrand, -2, -2 + c);

		//6.) Pass updated structures
		uAfterAtoms = newUAtoms;
		oAfterAtoms = newOAtoms;
	}
	else { //if there is not sea leve rise

		uAfter = uBefore;
		oAfter = oBefore;
		uAfterAtoms = uBeforeAtoms;
		oAfterAtoms = oBeforeAtoms;
	}
}



void disperse(const VecDoub &oBeforeAtoms, const VecDoub &uBeforeAtoms, const VecDoub &habitatQuality, const VecDoub &habitatQualityAtoms, 
	VecDoub &AtomsAfter, const int &indexOfSource, const double &b, const double &mu, const double &sigma, const double &rMax,
	const VecDoub &oBefore, const VecDoub &uBefore, VecDoub &oAfter, const VecDoub &densityDependenceVector){

    double modifier = .50;
	double locationOfSource = habitatQualityAtoms[indexOfSource];

	double densityDependence;

	if (uBeforeAtoms[indexOfSource] + oBeforeAtoms[indexOfSource] != 0){
		densityDependence = (uBeforeAtoms[indexOfSource]/(uBeforeAtoms[indexOfSource] + oBeforeAtoms[indexOfSource]));
	} else {
		densityDependence = 0;
	}

	double progeny = oBeforeAtoms[indexOfSource] + r(locationOfSource, mu, sigma, rMax)*oBeforeAtoms[indexOfSource]*densityDependence;
		
	int subintervals = habitatQuality.size();

    AtomsAfter[indexOfSource] += (1-modifier)*progeny; //remove A from the area in our atom if  it's going to be transferred
            // (Tails given by other atoms) + 90% of progeny
			// = tails from dispersal in the main growth function  - what was already there + 90 percent of progeny

	
	double lower, upper, lowerdiff, upperdiff, transferAmount;

	//first take care of the background
	for (int dest = 0; dest < subintervals; dest++){//go forwards to keep track of endpoints

		//subintervals are centered
		lower = habitatQuality[dest] - 0.5*(4.0/subintervals);
		upper = habitatQuality[dest] + 0.5*(4.0/subintervals);

		lowerdiff = lower - locationOfSource;
		upperdiff = upper - locationOfSource;

		transferAmount = (1/(2*b))*exp(-(abs(habitatQuality[dest]-locationOfSource)/b));

		oAfter[dest] += modifier*progeny*densityDependenceVector[dest]*transferAmount;
	
	} //end the for loop that goes destination by destination

    //Now, movement between atoms (only important atom is the one at -2, as an atom is a single point and won't receive area)

	if (uBeforeAtoms[0] + oBeforeAtoms[0] != 0){
		densityDependence = uBeforeAtoms[0]/(uBeforeAtoms[0] + oBeforeAtoms[0]);
	} else {
		densityDependence = 0;
	}
    
	AtomsAfter[0] += progeny*modifier*densityDependence*0.5*exp((-2.0 - habitatQualityAtoms[indexOfSource]) / b);
    //If the source is the atom at -2, then that atom is going to get up to half of modifier of progeny (depending on density dependence term)

    }

//Performs the growth step
//The size  of uAfterAtoms, oAfterAtoms and habitatQualityAtoms does not change
void grow(const VecDoub &uBefore, const VecDoub &oBefore, VecDoub &uAfter, VecDoub &oAfter, const double &mu,
	const double &sigma, const double &rMax, const double &b, const int &subintervals, const VecDoub &habitatQuality,
	const VecDoub &uBeforeAtoms, const VecDoub &oBeforeAtoms, VecDoub &uAfterAtoms, VecDoub &oAfterAtoms, const VecDoub &habitatQualityAtoms) {
		//1.) Reset  oAfter so that it is empty:
		for (int i = 0; i < subintervals; i++) {
			oAfter[i] = 0.0;
		}
		//Set AtomsAfter to 0 (only for o)
		for (int j = 0; j < habitatQualityAtoms.size(); j++){
			oAfterAtoms[j] = 0;
		}
		//oAfterAtoms = oBeforeAtoms;

		//Set up density dependent dispersal coeficcients 
		VecDoub densityDependenceVector(subintervals);
		setDensityDependenceVector(oBefore, uBefore, oBeforeAtoms, uBeforeAtoms, habitatQuality, habitatQualityAtoms, densityDependenceVector);

	 //2.) Create the integrand (this does not include the kernel)
	 double progeny, lower, upper, transferAmount, upperdiff, lowerdiff, progenyDensDep;

	//3.) Go one interval at a time, calculating the number of progeny in the interval
	for (int source = 0; source < subintervals; source++) {

		if (uBefore[source] + oBefore[source] != 0){
			progenyDensDep = (uBefore[source]/(uBefore[source]+oBefore[source]));
		} else {
			progenyDensDep = 0;
		}

        progeny = oBefore[source] + r(habitatQuality[source], mu, sigma, rMax)*oBefore[source]*progenyDensDep;

	 	//4.) Calculate how many of those progeny go to each of the k other intervals (dispersal kernel)
		//Use a density-dependent dispersal

	 	for (int dest = 0; dest < subintervals; dest++) {

            //subintervals are centered
		    lower = habitatQuality[dest] - 0.5*(4.0/subintervals);
		    upper = habitatQuality[dest] + 0.5*(4.0/subintervals);

            lowerdiff = lower - habitatQuality[source];
            upperdiff = upper - habitatQuality[source];

            transferAmount = 0.5*(signum(upperdiff)*(1 - exp(-abs(upperdiff) / b)) - signum(lowerdiff)*(1 - exp(-abs(lowerdiff) / b)));

	 		oAfter[dest] += progeny*densityDependenceVector[dest]*transferAmount;
	 		 //density-dependent dispersal

	 	}//End destination for loop

	 	//4.5) For each source, some goes to the atom at -2 and to the atom at 1 (with density-dependent dispersal)
        
        if (uBeforeAtoms[0] + oBeforeAtoms[0] != 0){
            double densityDependence;
		    densityDependence = (uBeforeAtoms[0] / (uBeforeAtoms[0] + oBeforeAtoms[0]));
            double correctionFactor = (4.0/subintervals); //Since we want to move probability, not probability density
	 	    oAfterAtoms[0] += (progeny*0.5*exp((-2.0 - habitatQuality[source]) / b) //receiving whole tail
	 		    * densityDependence)*correctionFactor; //density-depdent dispersal
        }
	// 	//No need to create an atom at two to get the right tail, because this atom would start with no unoccupied area, and would
	// 	//never gain any unoccupied area, because engineering can't occur past 1. The density-dependent dispersal would thereby
	// 	//prevent any occupied area from ever entering this atom, so may as well never create it in the first place

	 } //End source for loop


	// ////5.) Take care of all of the other atoms as the source, one at a time, using disperse left and disperse right
	//use density-dependent dispersal

	for (int asource = 0; asource < habitatQualityAtoms.size(); asource++){
        if (oBeforeAtoms[asource] + uBeforeAtoms[asource] != 0){
		    disperse(oBeforeAtoms, uBeforeAtoms, habitatQuality, habitatQualityAtoms, oAfterAtoms, asource, b, mu, sigma, 
		        rMax, oBefore, uBefore, oAfter, densityDependenceVector);
        }
	}

	 //6.) Accounting step

	finalize(oAfter);

	 for (int i = 0; i < subintervals; i++) {
        if (uBefore[i] + oBefore[i] < oAfter[i]){
	 	    oAfter[i] = uBefore[i] + oBefore[i];
        }
        uAfter[i] = uBefore[i] + oBefore[i] - oAfter[i];
	 }

	finalize(oAfterAtoms);
	//Kill everything in the atom at -2
	oAfterAtoms[0] = 0;
	  for (int i = 0; i < habitatQualityAtoms.size(); i++){
		  if (uBeforeAtoms[i] + oBeforeAtoms[i] < oAfterAtoms[i]){
			  oAfterAtoms[i] = uBeforeAtoms[i] + oBeforeAtoms[i];
		  }
	 	uAfterAtoms[i] = uBeforeAtoms[i] + oBeforeAtoms[i] - oAfterAtoms[i];

	  }
} //end of function


//A function describing the dependence of habitat modification on population density, used in engineering step
double h(const VecDoub &oBefore, const double &l, const double &sum, const int &subintervals, const double &hq, //habitat quality value
	const VecDoub &oBeforeAtoms, const VecDoub &habitatQualityAtoms, const VecDoub &habitatQuality) {
	if (sum <= 0.0){
		return 0.0;
	} else {
	//simpleIntegrand oIntegrand(habitatQuality, oBefore);
	double lower = MAX(hq - l, -2);
	double upper = MIN(hq + l, 2);
	//double areaOccupied = qromb(oIntegrand, lower, upper);
	double areaOccupied = 0;
	for (int i = 0; i < oBefore.size(); i++){
		if (habitatQuality[i] >= lower && habitatQuality[i] <= upper){
			areaOccupied += oBefore[i];
		}
	}

	areaOccupied *= (4.0/subintervals);

	//Add in atoms, if applicable
	for (int i = 0; i < habitatQualityAtoms.size(); i++){

		if (habitatQualityAtoms[i] <= upper && habitatQualityAtoms[i] >= lower){
			areaOccupied += oBeforeAtoms[i];
		}
	}

	return areaOccupied / sum;
	}
}

//Function describing the growth rate, which is needed for the population growth function, used in the growth step
double r(const double &x, const double &mu, const double &sigma, const double &rMax) {

	return rMax*(1 / pow(2 * PI*pow(sigma, 2.0), 0.5))*exp(-pow(x - mu, 2) / (2 * pow(sigma, 2.0)));

}

/*Function that performs interpolation*/
double interpolate(const VecDoub &X, const VecDoub &Y, const double x) {

	Poly_interp myfunc(X, Y, 4); //third order interpolation
	return myfunc.interp(x);
}

/*Converts NumericVector types to VecDoub types*/
VecDoub NumericVectorToVecDoub(const NumericVector &vec) {

	int sizeOfVec = vec.size();
	VecDoub outVec(sizeOfVec);
	for (int i = 0; i < sizeOfVec; i++) {
		outVec[i] = vec[i];
	}
	return outVec;
}

/*Converts VecDoub types to NumericVector types*/
NumericVector VecDoubToNumericVector(const VecDoub &vec) {

	int sizeOfVec = vec.size();
	NumericVector outVec(sizeOfVec);
	for (int i = 0; i < sizeOfVec; i++) {
		outVec[i] = vec[i];
	}
	return outVec;
}

/*Function that records each iteration to the NumericVector*/
void saveToRVector(const VecDoub &cStructure, NumericVector &rStructure, const int &pointer) {
	
		int end = cStructure.size();
	
		for (int i = 0; i < end; i++) {
			rStructure[pointer + i] = cStructure[i];
		}
	
	}


/*Fills in the time vector t*/
void setT(const int &subintervals, const int &timesteps, NumericVector &t){
	for (int ts = 0; ts <= timesteps; ts++) {
		for (int si = 0; si < subintervals; si++) {
			t[ts*subintervals + si] = ts;
		}
	}
}

//Sets up initial habitat quality, performed only once at beginning of program
void setHabitatQuality(VecDoub &habitatQuality, const int &subintervals) {
	for (int j = 0; j < subintervals; j++) {
		habitatQuality[j] = -2 + (j+0.5)*(4.0 / subintervals);
	}
}


//Set up x, performed only once at beginning of program
void setX(const int &subintervals, const int &timesteps, NumericVector &x, const VecDoub &habitatQuality) {
	for (int ts = 0; ts <= timesteps; ts++) {
		for (int si = 0; si < subintervals; si++) {
			x[ts*subintervals + si] = habitatQuality[si];
		}
	}
}


//get total number of Atoms for all of the timesteps added together
int getNumOfAtoms(const int &timesteps, const double &c) {
	
	int sum = 0;

	if (c == 0) {
		sum = 2 * (timesteps + 1);
	}
	else {
		int intervals = int(3.0 / c) - 1; //round down. Subtract 1 because we already have an atom at 1.

		// go timestep by timestep adding the number of atoms
		for (int i = 0; i <= timesteps; i++) {
			if (i < intervals) {
				sum += 2 + i; //at the initial timestep, there are already 2 atoms (at 1 and -2), and 1 more added each timestep
			}
			else { //once the interior region between -1 and 2 is filled with atoms, we just keep adding the same max amount of atoms every step
				sum += 2 + intervals;
			}
		}
	}
	return sum;
}

void setHabitatQualityAtoms(VecDoub &habitatQualityAtoms) {
	habitatQualityAtoms[0] = -2.0;
	habitatQualityAtoms[1] = 1.0;

}

void setUInitialAtoms(VecDoub &uInitialAtoms) {
	uInitialAtoms[0] = 0.0; //atom at -2
	uInitialAtoms[1] = 0.0; //atom at 1
}

void setOInitialAtoms(VecDoub &oInitialAtoms) {
	oInitialAtoms[0] = 0.0; //atom at -2
	oInitialAtoms[1] = 0.0; //atom at 1
}

void finalize(VecDoub &before){
	int subintervals = before.size();
	for (int i = 0; i < subintervals; i++){
		if (before[i] < 0){
			before[i] = 0;
		}
	}
}

//This is only done once at the beginning of the program
void setTAtoms(NumericVector &tAtoms, const int &timesteps, const double &c) {

	int numOfAtoms;
	int pointer = 0;
	// go timestep by timestep writing the timestep number the appropriate number of times
	for (int i = 0; i <= timesteps; i++) {
		if (c == 0) { //no sea level rise, always going to have the same number of atoms
			numOfAtoms = 2;
		}
		else { //if sea level rise, the number of atoms will increase
			int intervals = int(3.0 / c) - 1; //round down. Subtract 1 because we already have an atom at 1.
			if (i < intervals) {
				numOfAtoms = 2 + i; //at the initial timestep, there are already 2 atoms (at 1 and -2), and 1 more added each timestep
			}
			else { //once the interior region between -1 and 2 is filled with atoms, we just keep adding the same max amount of atoms every step
				numOfAtoms = 2 + intervals;
			}
		}
		for (int k = 0; k < numOfAtoms; k++) {
			tAtoms[pointer] = i;
			pointer++;
		}
	}
}

//This is only done once at the beginning of the program
void setXAtoms(NumericVector &xAtoms, const int &timesteps, const double &c) {
	int numOfAtoms;
	int pointer = 0;

	// go timestep by timestep writing the timestep number the appropriate number of times
	for (int i = 0; i <= timesteps; i++) { 

		//First, figure out how many atoms we have
		if (c == 0) {
			numOfAtoms = 2;
		}
		else {
			int intervals = int(3.0 / c) - 1; //round down. Subtract 1 because we already have an atom at 1.
			if (i < intervals) {
				numOfAtoms = 2 + i; //at the initial timestep, there are already 2 atoms (at 1 and -2), and 1 more added each timestep
			}
			else { //once the interior region between -1 and 2 is filled with atoms, we just keep adding the same max amount of atoms every step
				numOfAtoms = 2 + intervals;
			}
		}

		//Then fill in the appropriate values of habitat quality for each atom
		for (int k = 0; k < numOfAtoms; k++) {
			if (k == 0) {
				xAtoms[pointer] = -2.0; //the first habitat quality is -2
				pointer++;
			}
			else {
				xAtoms[pointer] = 1.0 - (numOfAtoms - k - 1)*c; 
				pointer++;
			}
		}
	}

}



double getTotalOccupiedArea(const VecDoub &oBefore, const VecDoub &oBeforeAtoms, const VecDoub &habitatQuality){

	double areaOccupied = 0;
	int subintervals = oBefore.size();
	for (int i = 0; i < subintervals; i++){
		areaOccupied += oBefore[i];
	}

	areaOccupied *= (4.0/subintervals);

	for (int i = 0; i < oBeforeAtoms.size(); i++){

		areaOccupied += oBeforeAtoms[i];
	}

	return areaOccupied;

}

double getDensityDependence(const VecDoub &oBefore, const VecDoub &uBefore, const VecDoub &oBeforeAtoms, const VecDoub &uBeforeAtoms,
const double &start, double &end, const VecDoub &habitatQuality, const VecDoub &habitatQualityAtoms){

	double areaOccupied = 0.0;
    double areaUnoccupied = 0.0;
    double subintervals = oBefore.size();

    for (int i = 0; i < subintervals; i++){
        if(habitatQuality[i] >= start && habitatQuality[i] <= end){
        areaOccupied += oBefore[i];
        areaUnoccupied += uBefore[i];
        }
    }
    areaOccupied *= (4.0/subintervals); //correction factor
    areaUnoccupied *= (4.0/subintervals);
		
        
        if (areaOccupied + areaUnoccupied == 0){
            return 0;
        } else {
		    return areaUnoccupied/(areaUnoccupied + areaOccupied);
        }
		
}


void setDensityDependenceVector(const VecDoub &oBefore, const VecDoub &uBefore, const VecDoub &oBeforeAtoms, const VecDoub &uBeforeAtoms, 
const VecDoub &habitatQuality, const VecDoub &habitatQualityAtoms, VecDoub &densityDependenceVector){

	int subintervals = habitatQuality.size();

	for(int dest = 0; dest < subintervals; dest++){
		if (oBefore[dest] + uBefore[dest] != 0){
			densityDependenceVector[dest] = uBefore[dest]/(oBefore[dest]+uBefore[dest]);
		} else {
			densityDependenceVector[dest] = 0;
		}

	}

}

/*Standard sign ("sgn") function that returns 1 if input is positive and -1 if input is negative*/
double signum(double &x) {
	return (x > 0) - (x < 0);
}




//Here is a function very similar to the main function, but this one only returns area occupied over time, not all of the data structures

// [[Rcpp::export()]]
List runModelForOccupied(const int &subintervals, const int &timesteps, const double &d, const double &m, const double &l,
	const double &mu, const double &sigma, const double &rMax, const double &b, const double &c, 
	const int &initialHabitatQuality, const bool &doEngineer, const bool &doErosion, const bool &doRiseSeaLevel,
	const bool &doGrow){
	
	//Find the number of atoms
	int numOfAtoms = getNumOfAtoms(timesteps, c); //number of atoms summed across all timesteps

	//Create C++ vectors

	//Create a vector to store values of habitat quality or "sediment level" and initialize this vector
	VecDoub habitatQuality(subintervals);
	setHabitatQuality(habitatQuality, subintervals);

	//Create a vector to store values of habitat quality for the finite number of atoms, and initialize this vector
	VecDoub habitatQualityAtoms(2); //2 is the number of initial atoms (one at -2 and one at 1)
	setHabitatQualityAtoms(habitatQualityAtoms);
	
	//Create vectors to store the initial distribution of unoccupied area and occupied area at each tidal height, and
	//initialize both of these vectors
	VecDoub uBefore(subintervals);
	setUInitial(uBefore, subintervals, initialHabitatQuality, habitatQuality);
	VecDoub oBefore(subintervals);
	setOInitial(oBefore, subintervals, habitatQuality);

	//Create vectors to store the initial vales of unoccupied and occupied area at each atom, and
	//initialize both of these vectors
	VecDoub uBeforeAtoms(2); //2 is the number of initial atoms (one at -2 and one at 1)
	setUInitialAtoms(uBeforeAtoms);
	VecDoub oBeforeAtoms(2); //2 is the number of initial atoms (one at -2 and one at 1)
	setOInitialAtoms(oBeforeAtoms);

	//After the growth step, we store the new area unoccupied and occupied at each tidal height in these vectors
	//There are 2 vectors for the distribution of values, and 2 vectors for the atoms
	VecDoub uAfter(subintervals);
	uAfter = uBefore; //In case none of the 4 steps are done, these "after" data structures still need to be set
	VecDoub oAfter(subintervals);
	oAfter = oBefore;
	VecDoub uAfterAtoms(2); //2 is the number of initial atoms (one at -2 and one at 1)
	uAfterAtoms = uBeforeAtoms;
	VecDoub oAfterAtoms(2); //2 is the number of initial atoms (one at -2 and one at 1)
	oAfterAtoms = oBeforeAtoms;

	//Declare R vectors
	NumericVector occupied(timesteps+1);

	//Find u and o, one timestep at a time:
	for (int i = 0; i < timesteps; i++) {

		cout << "Timestep " << i << endl; //Keep track of progress

		occupied[i] = getTotalOccupiedArea(oBefore, oBeforeAtoms, habitatQuality);
		if (i == 0){
			cout << "occupied: " << occupied[i] << endl;
		}

		//saveToRVector(habitatQuality, x, subintervals, i); Need to do this at the same time as the t vector for R

		//Perform all 4 steps
		if (doErosion){

			erode(uBefore, oBefore, uAfter, oAfter, d, subintervals, habitatQuality, uBeforeAtoms, oBeforeAtoms,
				uAfterAtoms, oAfterAtoms, habitatQualityAtoms);
			finalize(uAfter);
			finalize(oAfter);
			finalize(uAfterAtoms);
			finalize(oAfterAtoms);
			cout << "Erosion finished" << endl;

			//The vectors at the last step become the new vectors for the next timestep
			uBefore = uAfter;
			oBefore = oAfter;
			uBeforeAtoms = uAfterAtoms;
			oBeforeAtoms = oAfterAtoms;
		}
		if (doEngineer){

			engineer(uBefore, oBefore, uAfter, oAfter, m, l, subintervals, habitatQuality,
				uBeforeAtoms, oBeforeAtoms, uAfterAtoms, oAfterAtoms, habitatQualityAtoms);
			finalize(uAfter);
			finalize(oAfter);
			finalize(uAfterAtoms);
			finalize(oAfterAtoms);
			cout << "Engineering finished" << endl;

			//The vectors at the last step become the new vectors for the next timestep
			uBefore = uAfter;
			oBefore = oAfter;
			uBeforeAtoms = uAfterAtoms;
			oBeforeAtoms = oAfterAtoms;
		}
		if (doRiseSeaLevel){

			riseSeaLevel(uBefore, oBefore, uAfter, oAfter, c, subintervals, habitatQuality,  uBeforeAtoms, oBeforeAtoms,
				uAfterAtoms, oAfterAtoms, habitatQualityAtoms, i);
			finalize(uAfter);
			finalize(oAfter);
			finalize(uAfterAtoms);
			finalize(oAfterAtoms);
			
			cout << "Sea level rise finished" << endl;

			//The vectors at the last step become the new vectors for the next timestep
			uBefore = uAfter;
			oBefore = oAfter;
			uBeforeAtoms = uAfterAtoms;
			oBeforeAtoms = oAfterAtoms;
		}
		if (doGrow){

			grow(uBefore, oBefore, uAfter, oAfter, mu, sigma, rMax, b, subintervals, habitatQuality,
				uBeforeAtoms, oBeforeAtoms, uAfterAtoms, oAfterAtoms, habitatQualityAtoms);
			finalize(uAfter);
			finalize(oAfter);
			finalize(uAfterAtoms);
			finalize(oAfterAtoms);

			cout << "Growth finished" << endl;

			//The vectors at the last step become the new vectors for the next timestep
			uBefore = uAfter;
			oBefore = oAfter;
			uBeforeAtoms = uAfterAtoms;
			oBeforeAtoms = oAfterAtoms;
		}
	}

	//Save the last iteration (first the non-atoms, then the atoms)
	occupied[timesteps] = getTotalOccupiedArea(oAfter, oAfterAtoms, habitatQuality);

	//Return the data structure to R
	cout << "Returning the dataframe" << endl;
	return Rcpp::List::create(_["Occupied"] = occupied);
}
