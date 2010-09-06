
/*
 * mainHCAI_Com_IBM.cpp

 * Notes: just_means10a  fixes bugs with must_means10 that caused problems with the decolonization strategy (most importantly, numerical errors
 * were causing pdi, progprob and pdc values to change slowly over time).
 *
 * Not sure why it's called just_means, as the code should be good for outputting everything we need shouldn't it?
 * Julie? Any idea?
 *
 * May  18 2010
 * Authors: Julie Robotham & Ben Cooper
 *  Infection model with decolonisation and isolation interventions
 *
 */

#include <iostream>
#include "HCAI_COM_IBM.h"
#include <string>
#include <vector>
#include <map>
#include <list>
#include <math.h>
#include <fstream>
#include <istream>




/////////////////////CHANGED INCLUDE BIT ***************************************
////////////////////////////////////////////////////////////////////////////////

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//use on the mac
//#include </opt/local/include/gsl/gsl_rng.h>
//use on the mac
//#include </opt/local/include/gsl/gsl_randist.h>

//#include </usr/local/include/gsl/gsl_rng.h>
//#include </usr/local/include/gsl/gsl_randist.h>



//#include </opt/local/include/gsl/gsl_rng.h>
//#include </opt/local/include/gsl/gsl_randist.h>
//////////// /////////////////////////////////////////////////////////////////////



using namespace std;

bool verbose=true; //if true print more output

float prop_C_on_ad = 0.05; //proportion colonized on admission to ICU for non-high risk patients
float  prop_high_risk=0.0; // proportion of patients who are high risk
float prop_C_on_ad_high_risk = 0.0; //realtive prevalence in high risk group
int ISOCAP =50; //max isolation capacity (by default very high, so no limit on isolation capacity)
double effect_of_ISO=1;   //proportional reduction in transmission rate from a case when isolated
double effect_of_ISO_SD= 0.0; //0.622  //SD of proportional reduction in transmission rate from a case when isolated
double effect_of_secISO= 0.365;  // same but for the seoncary isolation measure used when primary isolation method capacity is exceeded
double effect_of_secISO_SD= 0; //0.622 // SD  for the seoncary isolation measure used when primary isolation method capacity is exceeded

double Pdc = 0.00013;//0.15;//low prev setting - calculated by trial and error to match literature scenarios//0.15 - high prev setting
double Pdi = 0.00000035;// 0.0035;//low prev setting - calculated by trial and error to match literature scenarios//0.035;- high prev setting  ///REALLY NEED TO LOOK AT THIS!!!
double ProgProb = 0.02; //from St Thomas' - 0.02 //0.045;//this will need to be changed if we are assuming all C->I movements are progressions
double IQ=1;//get_IQ(rng);//ignore heterogeniety for now
float CC_SENSITIVITY_MEAN=0.6815; //conventional culture
float CC_SENSITIVITY_SD=0;//0.194; //conventional culture
float CC_SPECIFICITY_MEAN = 0.882;
float CC_SPECIFICITY_SD = 0;//0.063 ;
float CA_SENSITIVITY_MEAN =0.8255; //chromagar after full incubation
float CA_SENSITIVITY_SD =0;//0.0427; //chromagar after full incubation
float CA_SPECIFICITY_MEAN=0.8305;
float CA_SPECIFICITY_SD =0;//0.1772;
float CA_EARLY_SENSITIVITY_MEAN=0.6217; //chromagar after full incubation
float CA_EARLY_SENSITIVITY_SD =0;//0.1249; //chromagar after full incubation
float CA_EARLY_SPECIFICITY_MEAN = 0.9713;  //chromager after short incubation
float CA_EARLY_SPECIFICITY_SD =0;// 0.0417;
float PCR_SENSITIVITY_MEAN=0.8840;  //pcr
float PCR_SENSITIVITY_SD = 0;//0.0510;  //pcr
float PCR_SPECIFICITY_MEAN=0.8380;
float PCR_SPECIFICITY_SD = 0;//0.0474;

int MEANCOLDURATION=100; //recovery should now be possible after discharge
int MEANINFDURATION=7;
int MEANHOMEDURATION=150;


double PropTreatmentSuccessful = 0.5855; ////from 2 studies (MUP +CHX, 1day and 1 week probabilitie of decol averaged - Rohr and HArbarth)//proportion for whom (total) decolonization happens
int length_of_treatment=5;//treatment lasts 5 days then stopped
double DecolEffectonPdc_MEAN =0.35;
double DecolEffectonPdc_SD = 0;//0.07;
double DecolEffectonPdi_MEAN =0.33;
double DecolEffectonPdi_SD =0;//0.21;
double DecolEffectonIQ_MEAN =1;
double DecolEffectonIQ_SD =0;
double DecolEffectonProgProb_MEAN =0;//0.69;
double DecolEffectonProgProb_SD =0;//0.18;

//costs (eventually these should be read from a file_

 float costofbedday = 1353;// thsi is currently set at the mean  - for comparison of startegies, for overlal CE with uncertainty need to use dist -  get_bd_cost(rng);//
 float dailyincremcostofsinglebediso = 53.89; // includes cost of contact precautions
 float dailyincremcostofIW = 67.28; // includes cost of contact precautions
 float dailycostofcontactprec = 15.15;//per patient value  - for a large general hospital (for a tertiary hospital would be 14.40)
 float costofswab = 4.25; //FOR LARGE GENERAL HOSPITAL (for a tertiary hosptial would be 4.04)
 float costofposCC = 8.78; //assuming that false postive screens have the same costs as true positives (similarly for false negs)
 float costofnegCC = 5.44;
 float costofposCA = 8.01;
 float costofnegCA = 4.67;
 float costofPCR = 20.87;//this is for IDI_MRSA **!!!!! REMEMEBER THAT THE SPEC AND SPEC NEED TO MATCH THE TEST YOU ARE COSTING!!!!!!!
 float costofinfetiontreatment =530.58; //from Rahul's treatment protocol
 float costofdecol  = 47.40;//
 float QoLinICU = 0.66;

int policynumber=0; //control policy index read in from command line
int simsperparameterset=100; //read in from command line
int numparametersetssampled=1;//

//extern float DecolEffectonIQ;
//extern float DecolEffectonPdc;
//extern float DecolEffectonProgProb;
extern double PropTreatmentSuccessful;


//double sensitivity= .75;
//double specificity = .86;

//lines below commented out as not using these as calculating additional length of stay differently
//long unsigned int LoS ;
//long unsigned int add_LoS;
//vector<unsigned long int> LoS_vector;
//vector<unsigned long int> add_LoS_vector;//additional length of stay vector

//these parameters stored as vectors which are initialised from file. These vectors specify a distribution of values
//and simulation runs sample from these values to choose parameters for a specific run
vector<double> sens_vector; //sensitivities (i.e. sensitivity of swab to detect MRSA carriage at any patient site)
vector<double> spec_vector; //specificities (i.e. specificity of swab to detect MRSA carriage at any patient site)

vector<float> Pdc_vector;
vector<float> Pdi_vector;
vector<float> Progprob_vector;

//these parameters hold probilities whic are specified for a number of days (but which are assumed to be fixed). So  element i of
// each vector lists probability of outcome for that type of patient for the ith day of stay
vector<float> dis_prob_SUS;   //discharge probability for a susceptible
vector<float> death_prob_SUS;  //death probability for a susceptible
vector<float> dis_prob_INF; //discharge probability for an infected
vector<float> death_prob_INF; //death probability for an infected

extern vector<float> home_duration_SUS; //duration at home for someone discharged susceptible (in practice every one's duration)

//vector<screeningpolicycomponent> weekdayscreening[7]; // screening policies for each day of week
//vector<screeningpolicycomponent> weeklypostadmission; // screening policies for screens every wk after a patient is admitted
//vector<screeningpolicycomponent> dischargescreening;


vector<unsigned long int> age_vector;
unsigned long int age;

//proba to probe are probabilities of being in one of 5 age groups
// (<19, 20-39, 40-59, 60-79, or 80+)
//note that they are cumulative probabilities and specify the probability
//if being in that age group or younger
//not using these at the moment, as not using age ------------
float proba; // probability  <19
float probb; // probability <40
float probc; // probability <59
float probd; // probability <79
float probe; // probability < infinite (i.e. 1)



double a; //total number of individuals in age group 0-19 in age distribution (from file)
double b; //total number of individuals in age group 20-39 in age distribution (from file)
double c; //total number of individuals in age group 40-59 age distribution (from file)
double d; //total number of individuals in age group 60-79 in age distribution (from file)
double e; //total number of individuals in age group 80+ in age distribution (from file)


//-----------------------------------------------------------


int numsim;



int main (int argc, char** argv) {
    //cout<<"mecamip ICU simulator: justmeans10 (May 2010)\n";
    gsl_rng * randgen; //pointer to gsl random number generator
    const gsl_rng_type * randgentype;

    char* Pdcfile=NULL;
    char* Pdifile=NULL;
    char* Progprobfile=NULL;


    randgentype =gsl_rng_mt19937;
    randgen = gsl_rng_alloc(randgentype);

    /*  *********** Process command line arguements ************    */
	//syntaix is: programname [optional arguments that start with "-" listed below] n1 n2 n3
	//where n1 is an integer specifying the policy number (policynumber)
	// n2 is number of simulations for each sampled set of parameter values (simsperparameterset)
	//n3 is number of sampled parameters (numparametersetssampled)


    while( (argc>1) && (argv[1][0]=='-')){

		switch(argv[1][1]){
			case 'v':
				// verbose  - if given print more output
				verbose=true;
				break;

			case 'p':
				//set prevalence on admission
				prop_C_on_ad=atof(&argv[1][2]);
				break;

			case 'r':
				//set relative prevalence in high risk group on admission
				prop_C_on_ad_high_risk=atof(&argv[1][2]);
				break;

			case 'd':
				//set duration of colonization in units of timestep
				MEANCOLDURATION= atoi(&argv[1][2]);
				break;

			case 'D':
				// set  in units of timestep
				MEANINFDURATION= atoi(&argv[1][2]);
				break;

			case 'h':
				//set proportion high risk
				prop_high_risk=atof(&argv[1][2]);

				break;

			case 'I':
				//set isoeffect mean
				effect_of_ISO=atof(&argv[1][2]);

				break;

			case 'i':
				//set isoeffect sd
				effect_of_ISO_SD=atof(&argv[1][2]);
				break;

			case 'J':
				//set secondary isoeffect mean
				effect_of_secISO=atof(&argv[1][2]);

				break;

			case 'j':
				//set secondary isoeffect sd
				effect_of_secISO_SD=atof(&argv[1][2]);
				break;

		//	case 'C':
				//set isocap (max isolation capacity)
			//	IsoCap=atoi(&argv[1][2]) ;
			//	break;

			case 'S':
				//set conventional culture sensitifity mean
	       		       	CC_SENSITIVITY_MEAN=atof(&argv[1][2]);
				break;

			case 's':
				//set conventional culture sensitifity sd
	     	       		CC_SENSITIVITY_SD=atof(&argv[1][2]);
				break;

			case 'F':
				//set conventional culature specificity mean
       			        CC_SPECIFICITY_MEAN=atof(&argv[1][2]);
				break;

			case 'f':
				// set conventional culature specificity sd
       				CC_SPECIFICITY_SD=atof(&argv[1][2]);
				break;

			case 'A':
				//set chromager (full) sensitifity mean
       				CA_SENSITIVITY_MEAN=atof(&argv[1][2]);
				break;

			case 'a':
				//set chromagar  culture sensitifity sd
       				CA_SENSITIVITY_SD=atof(&argv[1][2]);
				break;

			case 'E':
				//set chromager (full) specificity mean
       				CA_SPECIFICITY_MEAN=atof(&argv[1][2]);
				break;

			case 'e':
				//set chromagar  culture sensitifity sd
       				CA_SPECIFICITY_SD=atof(&argv[1][2]);
				break;

			case 'Y':
				//set chromager early sensitifity mean
      			        CA_EARLY_SENSITIVITY_MEAN=atof(&argv[1][2]);
				break;

			case 'y':
				//set chromagar early  culture sensitifity sd
				CA_EARLY_SENSITIVITY_SD=atof(&argv[1][2]);
				break;

			case 'Z':
				//set chromager early sensitifity mean
      			        CA_EARLY_SPECIFICITY_MEAN=atof(&argv[1][2]);
				break;

			case 'z':
				//set chromagar early  culture sensitifity sd
				CA_EARLY_SENSITIVITY_SD=atof(&argv[1][2]);
				break;

			case 'M':
				//set chromager early sensitifity mean
      			        PCR_SENSITIVITY_MEAN=atof(&argv[1][2]);
				break;

			case 'm':
				//set chromagar early  culture sensitifity sd
				PCR_SENSITIVITY_SD=atof(&argv[1][2]);
				break;

			case 'N':
				//set chromager early sensitifity mean
      			        PCR_SPECIFICITY_MEAN=atof(&argv[1][2]);
				break;

			case 'n':
				//set chromagar early  culture sensitifity sd
				PCR_SENSITIVITY_SD=atof(&argv[1][2]);
				break;


	     	       	case 'T':
				//set turnaround time (TAT) in units of timesteps
				cerr<<"Command line option T is not currently implemented" <<argv[1] << "\n";
				break;


			case 'w':
				//set delay for clinical swabs in units of timesteps
				//not currently implemented  - defaults to one time step
				cerr<<"Command line option w is not currently implemented" <<argv[1] << "\n";

				break;

			case 'g':
				//set number of negative swabs until considered clear
				//numbernegativeuntilconsideredcleared defaults to 3 swabs
				cerr<<"Command line option g is not currently implemented" <<argv[1] << "\n";


				break;

			case 'o':
				//set probability topical treatment clears carriage
				PropTreatmentSuccessful=atof(&argv[1][2]); //convert to float
				//may be better to use   std::stringstream s_number(&argv[1][2]);

				break;

			case 'B':
				//set decolonization effect on IQ
				DecolEffectonIQ_MEAN =atof(&argv[1][2]);
				break;

			case 'b':
				//set decolonization effect on IQ
				DecolEffectonIQ_SD =atof(&argv[1][2]);
				break;

			case 'K':
				//set decolonization effect on Pdc
				DecolEffectonPdc_MEAN =atof(&argv[1][2]);
				break;

			case 'k':
				//set decolonization effect on Pdc
				DecolEffectonPdc_SD =atof(&argv[1][2]);
				break;

			case 'U':
				//set decolonization effect on Pdc
				DecolEffectonPdi_MEAN =atof(&argv[1][2]);
				break;

			case 'u':
				//set decolonization effect on Pdc
				DecolEffectonPdi_SD =atof(&argv[1][2]);
				break;


			case 'X':
				//set decolonization effect on Progprob
				DecolEffectonProgProb_MEAN =atof(&argv[1][2]);
				break;

			case 'x':
				//set decolonization effect on Progprob
				DecolEffectonProgProb_SD =atof(&argv[1][2]);
				break;


			case 't':
				// set number of timesteps
				//not yet implemented...maybe better to do this with an env variable.
				// currently timesteps per day is set to 1 and defined by  DAILYTIMESTEPS
				cerr<<"Command line option t is not currently 	" <<argv[1] << "\n";
				break;

			case 'P':  //in this case we look at next 2 characters - if '1:' read in Pdc file, if '2:' read in Pdi; if '3:' read in ProgProbfile
				switch(argv[1][2]){
					case '1':
						if(argv[1][3]==':') {Pdcfile = &argv[1][4];}
						break;
					case '2':
						if(argv[1][3]==':') {Pdifile = &argv[1][4];};
						break;
					case '3':
						if(argv[1][3]==':') { Progprobfile = &argv[1][4];};
						break;
				}
				break;


			default:
				cerr<<"Unregonised command line option" <<argv[1] << "\n";

		}//end switch
		++argv; //move to next argument
		--argc; //one less argument to process



    } //end while
	// after reading in command line options should read in i) policy number to use ii) number of sims per parameter set iii) number of parameter sets
	if(argc!=4) {
		cerr<<"Uh oh. There's a problem: you need to supply 3 arguments (in addition to optional arguments): policy number; number of sims per parameter set; number of parameter sets \n";
		exit(1);
	} else {
	   	policynumber=atoi(argv[1]);
		simsperparameterset=atoi(argv[2]);
		numparametersetssampled=atoi(argv[3]);
	}
    if(verbose){ //print parameters
		cout<<"policynumber"<< policynumber<<"\n";
		cout<<"simsperparameterset"<<simsperparameterset<<"\n";
		cout<<"numparametersetssampled"<<numparametersetssampled <<"\n";

		if(Pdcfile!=NULL)  cout<<"**********************"<<Pdcfile<<'\n';
		if(Pdifile!=NULL)  cout<<"**********************"<<Pdifile<<'\n';
		if(Progprobfile!=NULL)  cout<<"**********************"<<Progprobfile<<'\n';
		cout<<"prevalence on admission " <<prop_C_on_ad << "\n";
		cout<<"prevalence on admission for high risk" <<prop_C_on_ad_high_risk << "\n";
		cout<<"duration of colonization (in units of timestep)" <<   MEANCOLDURATION  << "\n";
		cout<<"duration of infection (in units of timestep)" <<   MEANINFDURATION  << "\n";
		cout<<" proportion high risk  " << prop_high_risk  << "\n";

		cout<<"isolation effect mean  " << effect_of_ISO  << "\n";
		cout<<"isolation effect sd   " << effect_of_ISO_SD  << "\n";
		cout<<" secondary isoeffect mean  " << effect_of_secISO  << "\n";
		cout<<" secondary isoeffect SD " << effect_of_secISO_SD  << "\n";
		cout<<" isolation capacity  " << ISOCAP  << "\n";
		cout<<"probability topical treatment clears carriage   " << PropTreatmentSuccessful  << "\n";
		//cout<<" decolonization effect on IQ  " << DecolEffectonIQ  << "\n";
		//cout<<" decolonization effect on Pdc   " <<  DecolEffectonPdc << "\n";
		//cout<<"  decolonization effect on Progprob  " <<  DecolEffectonProgProb   << "\n";
		//    cout<<"   " <<   << "\n";

    }

    /************ Get disrtibutions from files  ************    */


    read_in_files(Pdcfile, Pdifile, Progprobfile);


    //modify code to make sure the values of these read in get used.

	/**  -----titles at top of output file -------*/



    if(BATCHMODE>0){
		cout<< "susbedday "<<" Colbedday "<<" Infbedday "<<" Isobedday "<<" Decolcount "<<" no_ad_screens "<<" no_wkly_screens "<<" no_clin_screens "<<" no_pos_screens "<<" no_neg_screens ";
		cout<<" no_pos_CC_screens "<< " no_neg_CC_screens "<<" no_pos_CA_screens "<< " no_neg_CA_screens "<<" no_pos_CA_early_screens "<< " no_neg_CA_early_screens "<<"totalcoldischarged"<< " no_PCR_screens ";
		cout<<" cumulative_StoC "<<" cumulative_StoI ";
		cout<<" cumulative_CtoI " <<" cumulative_admissions "<<" cumulative_discharges "<< " cumulative_readmissions "<< " cumulative_deaths ";
		cout<< " cumulative_appisodays " <<" cumulative_inappisodays "<<" cumulative_unisodays " <<" cumulative_good_bd ";
		cout<<" beddaycosts "<<" isolationcosts "<<" decolcosts "<<" swabbingcosts "<<" screeningcosts "<<" treatmentcosts "<<" totalcosts ";
		cout<<" healthbenefitsinICU " <<"\n";//"healthbenefitspostdischarge"<<
    }

   long  int num_Pdc_values=Pdc_vector.size();
   long  int num_Pdi_values=Pdi_vector.size();
   long  int num_Progprob_values=Progprob_vector.size();


    for ( numsim=0; numsim<numparametersetssampled; numsim++){






      //if Pdc values were sampled, then choose one randomly
      if(num_Pdc_values>0) {
        int sample=int(gsl_rng_uniform(randgen)*num_Pdc_values);
        Pdc=Pdc_vector[sample];
        if(verbose)       cout<<" Sampled Pdc value: "<<Pdc<<"\n";
      }

      //if Pdi values were sampled, then choose one randomly
      if(num_Pdi_values>0) {
        int sample=int(gsl_rng_uniform(randgen)*num_Pdi_values);
        Pdi=Pdi_vector[sample];
        if(verbose) cout<<" Sampled Pdi value: "<<Pdi <<"\n";
      }


     //if Progprob values were sampled, then choose one randomly
      if(num_Progprob_values>0) {
        int sample=int(gsl_rng_uniform(randgen)*num_Progprob_values);
        ProgProb=Progprob_vector[sample];
        if(verbose) cout<<" Sampled Pprogprob value: "<<ProgProb<<"\n";
      }

      //create a new patients objects numparametersetssampled times (with randomly specified control effectivneess etc
      // based on sampled values - where distributions mean and SD can be sepcified o nthe commnaand samples paremter from thise


		patients *testpatient = new patients(randgen);//randgen

		testpatient->setpolicy(randgen,policynumber);

     //for each such parameter set perform   simsperparameterse iterations and output the mean
		testpatient->run_n_simulations(4380, 0, randgen, simsperparameterset);//(1825, 0, randgen);//stoptime, no initial cases WHICH SHOULD BE 0 AS WE ARE ADMITTING THEM INSTEAD, random no

		//   cout<<"ICU data after the simulation \n";
		//testpatient->print_ICU_data(); //(std::vector<current_ICU_patients> ICU_patients);

        testpatient->print_events();
		// testpatient->print_all_patients_data();

		//cout<<" Day "<<" SinICU "<<" CinICU "<<" IinICU "<<" ISOinICU "<<"/n ";




		delete testpatient;



      //  testpatient.print_all_data();
        //testpatient.perform_events(randgen);
    } //for numsim


	return 0;
}

