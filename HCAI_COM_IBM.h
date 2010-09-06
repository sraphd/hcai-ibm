

/*
 * just_means10a.h

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



#define DEFAULT_DEBUG 103
// if BATCHM0DE is set to greater than 0 BATCH MODE output is used...i.e. only one line of ouput per simulation and no daily info
#define BATCHMODE 1


// OUTPUTLEVEL specified how much info to input. 0 is minimum, 1 is more etc

#define OUTPUTLEVEL 0


//DEBUG LEVEL 1 = admissions and discharges
//DEBUG LEVEL 2 = screening
//DEBUG LEVEL 3 = admissions and discharges
//DEBUG LEVEL 4 = d screening
//DEBUG LEVEL 6 = discharge screening
//DEBUG LEVEL 8 = isolation
//DEBUG LEVEL 16 = prevalence output
//DEBUG LEVEL 32 = infection
//DEBUG LEVEL 64 = decolonisation
//DEBUG LEVEL 128 = read-in files checks
//DEBUG LEVEL 72 = readmission check
//DEBUG LEVEL 75 = home duration
//DEBUG LEVEL 100=follow one patient number
//DEBUG LEVEL 103=print out simulation to txt file


//all these defines need to be changed to input that can be changed at runtime
#define DEFAULT_POPSIZE 180000
#define DEFAULT_ICU_SIZE 1000
//#define START_TIME_OF_INT 00
#define DAILYTIMESTEPS 1

//#define stoptime 1825
// should eventually specify number of daily time steps in an environment variable

// #define STEPSIZE 1.0
//#define prop_scr 1.0  //set to x where 0<=x<=1. Less than 1 means some screens missing (eg.e .9 would meet 90% compliance).
//#define choose_decol 8.0

//#define sensitivity 1.0
//#define specificity 1.0

//IsoCap is isolation capacity -if more than this number are in isolation some have to wait
//needs to be changed to input that can be specified at run-time

//#define IsoCap 60

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <math.h>
#include <fstream>
#include <istream>


/////////////////////CHANGED INCLUDE BIT ***************************************
///////////////////////////////////////////////////////////////////////////////

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//#include </opt/local/include/gsl/gsl_rng.h>
//#include </opt/local/include/gsl/gsl_randist.h>


//#include </opt/local/include/gsl/gsl_rng.h>
//#include </opt/local/include/gsl/gsl_randist.h>

//#include </usr/local/include/gsl/gsl_rng.h>
//#include </usr/local/include/gsl/gsl_randist.h>



///////////////////////////////////////////////////////////////////////////////



extern  float costofbedday;
extern float dailyincremcostofsinglebediso;
extern float dailyincremcostofIW;
extern float dailycostofcontactprec;
extern float costofswab;
extern float costofposCC;
extern float costofnegCC;
extern float costofposCA;
extern float costofnegCA;
extern float costofPCR;
extern float costofinfetiontreatment;
extern float costofdecol;
extern float QoLinICU;


extern int numsim;


using namespace std;

extern gsl_rng * randgen; //pointer to gsl random number generator
extern	const gsl_rng_type * randgentype;

// Change line below to output different amounts of debuggin information


//*********************SPECIALTIES**********************************************
const char GM = 201;
const char ICU = 202;
const char COM = 203;
const char NA =204;

typedef char T_Specialty;


//*********************DISEASE STATES*******************************************
// types for disease state - store in one byte
//space saving excercise

const char SUSCEPTIBLE = 1;
const char COLONIZED = 2;
const char INFECTED =3;
const char COLANDINF =4;  //colonized and infected

//event types
const char COLONIZE = 5;
const char INFECT =6;
//const char INFECTFROMCOL = 103;
const char ICU_DIS_AND_ADMIT =7;
const char SCHEDULED_SCREEN = 8;
const char CLINICAL_SCREEN = 26;
const char FINISH_TREATMENT = 9;
const char ACT_ON_POS_SCREEN = 10;
//const char ACT_ON_RESCREEN = 111;
const char ACT_ON_NEG_SCREEN = 11;
const char DEATH = 12;
const char RECOVERFROMINF =13;
const char RECOVERFROMCOL =14;
const char READMISSION = 15;

//const char PROGRESSION = 106;

typedef char T_DisState ;	//disease states


//*********************DISCHARGE STATES*******************************************
// types for discharge state - store in one byte
//space saving excercise

const char ALIVE = 15;
const char DEAD = 16;


typedef char T_Discharge ;	//dicharge states

//*********************AWARENESS STATES*******************************************

const char UNKNOWN = 17;
const char KNOWN_S = 18;
const char KNOWN_C = 19;
const char KNOWN_I = 20;
const char BELIEVED_S_ACTUALLY_C = 21;
const char BELIEVED_S_ACTUALLY_I = 22;
const char BELIEVED_C_ACTUALLY_S = 23;

typedef char T_AwarenessState;

//*********************SCREEN RESULTS*******************************************

const char POSITIVE = 1;
const char NEGATIVE = 2;
const char NORESULT = 3;

//*********************RISK GROUPS**********************************************
const char HIGH = 24;
const char LOW = 25;

typedef char T_RiskGroup;

//**********************INTERVENTIONS*******************************************
//types for interventions - also stored as characters to save space.
//Can have combinations of these



const char NONE=0;
const char CC= 32 ; //conventional culture for all
const char CA = 33; //chromagar - with full incubation
const char CA_early =34;//chromagar - at early incubation time
const char PCR =35;
const char I =36; // hypothetical ideal test



typedef char screeningtype; //should take one of the above tyes



/** *********************declarations of KEY SIMULATION PARAMETERS*******************************************/

// All thise stuff moved so these are now defined in main


//default values defined here...but need to be able to overwrite these by eeither reading from a file, specifying on the command line etc

    const short int timestepsperday=DAILYTIMESTEPS; //defaults to daily timesteps...but we need to allow for other possibilities
                                       //note that all parameters such as TAT must be expressed in units of timesteps
                                       //for patient admission and discharge events - which are only estimated on a daily basis
                                       //we could just specify that admissions and discharges only occur on the first time step of each day.
//double sensitivity;
//double specificity;

    /** transmission parameters*/
// @@@ modify these so pdc pdi and progprob are are read in from a file
//  double xxxxxxxxxxtestxxxxxxxxxxxx;
//   double Pdc = 0.016;//0.15;//low prev setting - calculated by trial and error to match literature scenarios//0.15 - high prev setting
//    double Pdi = 0.001;// 0.0035;//low prev setting - calculated by trial and error to match literature scenarios//0.035;- high prev setting
//    double ProgProb = 0.18; //from St Thomas' - 0.02 //0.045;//this will need to be changed if we are assuming all C->I movements are progressions
//    const double IQ=1;//get_IQ(rng);//ignore heterogeniety for now
//    int MEANCOLDURATION=370; // duration of colonization in units of timesteps (should be possible to specify this at run time)
//    int  MEANINFDURATION=370; //30; // mean duration of infection in units of timesteps(should be possible to specify this at run time)

    /** intervention parameters - screening*/
    //
    /** intervention parameters - isolation*/

    //const double ISOEFFECT;//=1.0; //0.35;// this should be constant in each simulation, but drawn from a distribution
    //const double SECISOEFFECT;//=0.1; /////COMPLETELY MADE UP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    /** intervention parameters - decolonization*/

   /*  double PropTreatmentSuccessful = 0.5855; ////from 2 studies (MUP +CHX, 1day and 1 week probabilitie of decol averaged - Rohr and HArbarth)//proportion for whom (total) decolonization happens */
/*     int length_of_treatment=5;//treatment lasts 5 days then stopped */
/*     float DecolEffectonPdc = 0.52;//from Ridenour et al - BUT this is Mup + CHX// 0.72 Barry's guess */
/*     float DecolEffectonIQ =0;// assume that all of the effect on IQ is captured by reversion to S after treatment is successful //0.72;//from Barry's guess */
/*     float DecolEffectonProgProb =0.9;//from Muller et al. */






/**STRUCTURES*/

struct all_patients{
	unsigned long int patientid; //individual patient identification number
	//(NB: objects of type unsigned are incapable of containing negative values)
	unsigned long int age;
	char gender;
    short int ICU_no;
	short int no_days_in_hos;
	short int no_days_at_home;
	short int no_hos;
	bool everpositiveswabthisadmission;   //true if any positive swab so far this admission
	unsigned short int timecolonized; //timestep colonized, if ever colonized
	unsigned short int timeinfected; //timestep infected, if ever infected
    unsigned short int timerecovered; //timestep recovered, if ever recovered
	double Pdc, Pdi;

    //Pdc = actual daily risk of transmission leading to colonization for that individual if there is 1 infectious patient
	//Pdi = actual daily risk of transmission leading to infection for that individual if there is 1 infectious patient

  	double ProgProb; //probablity of progressing from C to I
	double IQ; //infectivity quotient  - how infectious somebody is between 0 and 1. Chnges with intervention - at the moment the same for C or I

	double baselinePdc, baselinePdi, baselineProgProb, baselineIQ;
	//baselinePdc,  baselinePdi, baselineProgProb and baselineIQ are there to store baseine values without interventions, so that when interventinos Pdc and Pdi can
	// return to these values.

    bool isolation; //true if patient is currently isolated (using whatever isolation means are currently being ued)
    bool secndisolation; // true if patient should be isoalted but the isoaltion facility is at capacity then put under secondary isolation e.g. cohorting
    bool decolonization; //true if patient is undergoing treatment (true for the whole treatment period)
    bool readmission; //true if patient is scheduled to be readmitted

   // bool decolonization_mup; //decolonization with mupirocin currently being used
   // bool decolonization_mup_chx; //decolonization with mupirocin & chlorhexidine currently being used
    bool nursecohorting; // nurse cohorting currently being used for this patient

//	T_Intervention interventions; // interventions currently in place, as defined in the 'patients' class
//  only reason for storing multiple interventions in the same varaible would be to save space. That was done in the flu
//  household model because we wanted to model several million individuals and memory was a limitation
//  in this case we have only a few thousand and memory is not a limiation so can use a separate variable for each component of an intervention
//  currently in place
    char previous_screen_result; //POSITIVE, NEGATIVE or NA
    short int num_consec_neg_screens_following_a_pos; //used to evaluate whether patient should be considered free of colonizatino/infection
	T_DisState disease_state; // disease state is defined in 'patients' class
    T_AwarenessState awareness_state; //what the HCAI status of the patient is believed to be following screening
    T_Specialty specialty;
    T_Discharge discharge_state;//whether the individual is discharged alive or dead
    T_RiskGroup risk_group;//i.e. high or low risk of colonization on admission


    /**
    //not using this now that discharges are calculated daily
    unsigned short int extra_stay_to_be_added; //flag to check (on discharge) whether the individual has become col or inf - and if so, add on some extra stay
    //1 = became col from sus
    //2 = became inf from sus
    //3 = became inf from col
    */
};

struct screening_parameters{
 double sensitivity; //sensitivity of screen for detetcitng MRSA at any site (betwen 0 and 1)
 double specificity; //specificity of screen for detetcitng MRSA at any site (betwen 0 and 1)
 unsigned short int tat; //turnaround time from screening (in units of time steps) from screen being taken to results acted on
};

struct screeningpolicycomponent{
 screeningtype technology;   // type of of culture (and possibly swab) used for screening
 bool targeted;  //if true for screening applies only to high risk patients
 bool previouslypositive; //if true policy applies only to patients who were reviously positive
};

struct interventionpolicy{
  // data structure to define intervention and screening policy
  vector<screeningpolicycomponent> admissionscreening;
  vector<screeningpolicycomponent> weekdayscreening[7]; // screening policies for each day of week
  vector<screeningpolicycomponent> weeklypostadmission; // screening policies for screens every wk after a patient is admitted
  vector<screeningpolicycomponent> dischargescreening;
  vector<screeningpolicycomponent> clinicalscreening;
  //type of admission screening (could be none, conv culture for all, chromagar for all,
  //pcr for all, chromagar with two readings (one early), targeted versions of the aforementioned etc
  // or combinations of any of the above (hence stored as a vector - each elelement is one comnent

  bool screeninghighriskonly;

  bool preemptiveisolationforall; //isolate all patients on admission until negative
  bool preemptiveisolationforhighrisk;//isolate all patients on admission until confrimed negative

  bool blanketdecolforall; //put all patients on decolonisation regime on admission
  bool decolforhighrisk;
  bool decolifpositive;


  unsigned short int numbernegativeuntilconsideredcleared;//number of consecutive -ve swabs required to assume clear of MRSA after being postive
  bool isolateifpositive; //isolate only MRSA positives
  int timesteptoimplement; //at what time step does the policy kick in
  double primary_isolation_effectiveness;

  int primary_isolation_capacity;


  double secondary_isolation_effectiveness;
  double bodywash_pdc_effectiveness;
  double bodywash_pdi_effectiveness;
  double bodywash_IQ_effectiveness;
  double bodywash_progprob_effectiveness;
  float proportionscreenedonadmission; //proportion of target population  screened on admission
  float proportionscreenedweekly; //proportion of target population  screened on a weekly or scheduled screen
  float proportionscreenedondischarge; //proportion of target population  screened on discharge
  int delaybeforeclinicalswab; //in timesteps

};


struct current_ICU_patients{

    all_patients* ICU_patient_ptr; //pointer to first patient in ICU
    unsigned short int bed_no;

	//bool in_sideroom; //in sideroom or not

};


struct events{  //event to be carried out some time in the future

    unsigned char event_type; //code indicating which event to carry out (defined in 'patients' class)
    all_patients* patient_ptr; //pointer to one patient

};


/*struct durations{  //event to be carried out some time in the future

    unsigned char duration_length; //code indicating which event to carry out (defined in 'patients' class)
    all_patients* patient_ptr; //pointer to one patient

};
*/


/*struct ICUstate {

    int time;
    long int SinICU;
    long int CinICU;
    long int IinICU;
    long int ISOinICU;
    long int noICU;
    long int HIGHRISKinICU;

};*/




void read_in_parameters (gsl_rng * rng);
void read_in_files (char*, char*, char*);
void pick_patient_function();




//double effect_of_ISO (gsl_rng * rng);
//double effect_of_secISO (gsl_rng * rng);
double get_IQ (gsl_rng * rng);
double get_bd_cost (gsl_rng * rng);
double get_CC_sensitivity (gsl_rng * rng);
double get_CC_specificity (gsl_rng * rng);
double get_CA_sensitivity (gsl_rng * rng);
double get_CA_specificity (gsl_rng * rng);
double get_CA_early_specificity (gsl_rng * rng);
double get_CA_early_sensitivity (gsl_rng * rng);
double get_PCR_specificity (gsl_rng * rng);
double get_PCR_sensitivity (gsl_rng * rng);
//float get_primary_isolation_effectiveness(gsl_rng * rng);
//float get_secondaryary_isolation_effectiveness(gsl_rng * rng);

int get_isocap(void);
double  get_effect_of_secISO(gsl_rng * rng);
double  get_effect_of_ISO(gsl_rng * rng);
double get_pdc_effect_of_bodywash(gsl_rng * rng);
double get_pdi_effect_of_bodywash(gsl_rng * rng);
double get_IQ_effect_of_bodywash(gsl_rng * rng);
double get_progprob_effect_of_bodywash(gsl_rng * rng);


/** PATIENTS CLASS*/

class patients{

public:

   ~patients(void);
	patients(gsl_rng* rng);//constructor of patients that takes a random number
	//~patients(); //DESTRUCTOR
	unsigned long int get_number_susceptible();
	unsigned long int get_number_colonized();
	unsigned long int get_number_infected();
    unsigned long int get_number_colonized_and_infected();

	void reset_all_counts(); //resets all the outcome variables to zero, enabling a new simulation run

	void print_all_data();//declaration of print functions
	void print_ICU_data();

	void process_ICU_admissions(gsl_rng * rng);//movement functions
    void process_ICU_discharges(all_patients* patient_ptr,gsl_rng *rng );
    void process_ICU_deaths(all_patients* patient_ptr, gsl_rng * rng);
    void process_readmission(gsl_rng * rng);

    void admission_screening(gsl_rng *rng, int pick_patient);//intervention functions
    void discharge_screening(all_patients* patient_ptr, gsl_rng *rng);//intervention functions
    void scheduled_screening( all_patients* patient_ptr, gsl_rng *rng);
    void clinical_screening( all_patients* patient_ptr, gsl_rng *rng);
    //above procudure performs scheduled screening and schedules a new screen for the same patient at a time screeninginterval in the future

    void scheduled_dayofweek_screening( all_patients* patient_ptr, gsl_rng *rng);
    //above procedure performs a scheduled screen on a specific day of the week and doesn't schedule any new screens

    void implement_control_measures_for_positives(int pick_patient, gsl_rng *rng); //apply whatever control measures are in place
    void implement_control_measures_for_positives(all_patients* patient_ptr, gsl_rng *rng);//apply whatever control measures are in place
    void remove_control_measures_for_those_believed_negative(all_patients* patient_ptr, gsl_rng *rng);
    void implement_control_measures_for_all(gsl_rng *rng, int pick_patient);
    void implement_control_measures_for_high_risk(gsl_rng *rng, int pick_patient);
    //void implement_control_measures_for_all(all_patients* patient_ptr, gsl_rng *rng);//apply whatever control measures are in place

    void remove_control_measures_for_high_risk(all_patients* patient_ptr, gsl_rng *rng);

    void decol(all_patients* patient_ptr, gsl_rng *rng);
    void decol(gsl_rng *rng, int pick_patient);
    //void CHXdecol(all_patients* patient_ptr, gsl_rng *rng);
    //void decolonization(gsl_rng *rng, int pick_patient);
    //void decolonization(all_patients* patient_ptr, gsl_rng *rng);
    void end_decolonization(all_patients* patient_ptr, gsl_rng *rng);
    void end_decolonization(gsl_rng *rng, int pick_patient);


    void perform_events(gsl_rng *rng); //performs all scheduled events at the end of a time step
    void print_events(); //performs all scheduled events at the end of a time step

    void app_or_inapp_isolation_days();
    void bed_day_counts();//counting routines - this updates cumulative values of coloinzed, infected, susceptible and isolated bed days for the current simulation
    void LOS_counts(); //just keeps trackof number time of days each patient has spend in the ICU so far;
    void LOH_counts(); //just keeps trackof number time of days each patient has spend at home after discharge so far;



    std::vector<unsigned long int> iso_queue; //really this shoudl be initliased
    std::vector<unsigned long int> admis_queue;



    int pick_LoS(gsl_rng* rng);
    int pick_add_LoS(gsl_rng* rng);
    int calculate_characteristics(gsl_rng* rng);


	T_DisState get_disease_state(const long int patientid);//{ //gets person's disease state
//	  return(hos_pop[patientid].disease_state);
//	}


//	void set_intervention(const long int patientid, T_Intervention action){//applies control measure to person
//	  hos_pop[patientid].interventions=action|hos_pop[patientid].interventions; //inclusive OR...so adds interventions to those that exist
//	  //The bitwise inclusive OR operator (|)
//	  //compares each bit of its first operand to the corresponding bit of its second operand.
//	}
//
//
//	T_Intervention get_intervention(const long int patientid){//get control measure applied to person
//	  return(hos_pop[patientid].interventions);
//	}
//


	void run_n_simulations(const int stoptime, const int numinitialcases, gsl_rng *rng, const long int n);

	void print_all_patients_data();

	void setpolicy(gsl_rng* rng, const int); //set interventionpolicy & screening parameter etc
        // attempt to do this first by looking at the command line options -
        //if nothing there look to environment variables and if nothing there ask user to enter

    void printpolicy(); // print current intervention policy

    float getbeddaycosts();
	float getisolationcosts();
	float getdecolcosts();
	float getswabbingcosts();
	float getscreeningcosts();
	float gettreatmentcosts();

	float gethealthbenefitsinICU();
     //double getsensitivity(screeningtype); //should return senstivity of  the screening type

     //double getspecificity(screeningtype); //should return specificity of  the screening type

    // double gettat(screeningtype); //should return tat of  the screening type

     bool isthereweeklypostadmissionscreening(void);//true if there is, false if there isn't

//     bool isthereweeklypostadmissionscreening(void){ //true if there is, false if there isn't
//       if(policy.weeklypostadmission.begin()==policy.weeklypostadmission.end()){
//         return(false); //since no weekly post admsisin screening
//       } else {
//         return(true);
//       }
//    };

    int get_col_duration(gsl_rng * rng); //returns duration of colonization in units of timestes- selected randomly from an exponential distribution
    int get_home_duration(gsl_rng * rng); //returns duration of home stay between readmissions
    int get_inf_duration(gsl_rng * rng); //return returns duration of infection in units of timestes- selected randomly from an exponential distribution
    int get_weekday_number(void);//returns 0 for monday, 1 for tuesday etc, where time=0 corresponds to Monday (accounts for timsteps per day (assumed to be >=1)
    bool is_it_first_time_step_of_day(void);// returns true if yes false if not
private:


    unsigned short int debug; //paramter to set de-bugging level, defaults to DEFAULT_DEBUG

	unsigned short int time; //
	const int  timestepsperday; // number of timesteps in one day = 1/timestep size in days. Must be >=1
	const unsigned long int popsize;  //number of people in population, defaults to  DEFAULT_POPSIZE

	struct all_patients hos_pop[DEFAULT_POPSIZE];
	// 'hos_pop' is the vector of all patients who ever go into ICU
	//of struct 'all_patients'

	struct current_ICU_patients ICU_patients[DEFAULT_ICU_SIZE];

    std::map<unsigned short int, std::vector<events> > future_movement_events; //pending admission/discharge/death events
    std::map<unsigned short int, std::vector<events> > future_infection_events; //pending infection events
    std::map<unsigned short int, std::vector<events> > future_intervention_events; //pending events related to intervention measures

    std::map<int, unsigned long int> duration_list; //list of home duration assigned
    //std::map<unsigned short int, std::vector<events> > future_screening_events;

    //intervention parameters
    struct screening_parameters conventional_culture; //holds parameters for conventional culture
    struct screening_parameters chromagar; //holds parameters for chromagar (after full incubation)
    struct screening_parameters chromagar_early; //holds parameters for early reading from chromagar
    struct screening_parameters pcr; //holds parameters for pcr screening
    struct screening_parameters ideal; //holds parameters for pcr screening
    struct interventionpolicy policy; //defintion of the screening and intervention policy

    float primary_isolation_effectiveness;
    float secondary_isolation_effectiveness;

	//summary variables
	long int S; //number susceptible
	long int C; //number colonized
    long int I; //number infected
	long int CI; //number colonized and infected
	long int ISO; //number in isolation
	long int SECISO;//number in secondary isoaltion - when primary isoaltion at capacity
	long int D; //number decontaminated ?????
	long int ComC; //number discharged colonised
	long int ReC;


    //ICU summary variables
 //  struct ICUstate;


   //ICU summary variables
   long int SinICU;
   long int CinICU;
   long int IinICU;
   long int ISOinICU;
   long int noICU;
   long int HIGHRISKinICU;

   long int numbersusceptible;
   long int numbercolonised;
   long int numberinfected;

//   ICUstate** stepresults;






    //bed day counts
    long int susbedday; //total number of susceptible bed days
    long int colbedday; //total number of colonized bed days
    long int infbedday; //total number of infected bed days
    long int isobedday; //total number of isolated bed days

    long int decolcount; //total number of decolonised bed days


    //isolation bed day counts
    long int appisodays; //cumulative number of appropraite isolation days
    long int inappisodays; //cumulative number of inappropriate isolation days
    long int unisodays; //cumulative number of unisolated days

    // will have to redo additional days due to infection
    long int addaysduetoinf;

    long int no_ad_screens;
    long int no_wkly_screens;
    long int no_clin_screens;
    long int no_pos_screens;
    long int no_neg_screens;

    long int no_pos_CC_screens;
    long int no_neg_CC_screens;
    long int no_pos_CA_screens;
    long int no_neg_CA_screens;
    long int no_pos_CA_early_screens;
    long int no_neg_CA_early_screens;
    long int no_PCR_screens;
    long int no_ideal_screens;


    long int dis_alive;
    long int dis_dead;

    long int cumulative_StoC; //total mumber of Susceptible to Colonized events in simulation so far
    long int cumulative_StoI; //total mumber of Susceptible to Infected events in simulation so far
    long int cumulative_CtoI; //total mumber of Colonized to Infected events in simulation so far
    long int cumulative_Conadmission;
    long int cumulative_admissions; //total mumber of admissions  in simulation so far
    long int cumulative_discharges; //total mumber of discharges  in simulation so far
    long int cumulative_deaths; //total mumber of deaths  in simulation so far
    long int cumulative_appisodays; //total mumber of appropriate isolation days  in simulation so far
    long int cumulative_inappisodays; //total mumber of inappropriate isolation days  in simulation so far
    long int cumulative_unisodays; //total mumber of inappropriate isolation days  in simulation so far
    long int cumulative_good_bd;//total  spent in hospital by those who are discharged alive
    long int totalcoldischarged; //total colonised on discharge
    long int cumulative_readmissions; //total number of readmissions in simulation so far


    //COSTS////////////////////////////////////////////


    float beddaycosts;
    float isolationcosts;
    float decolcosts;
    float swabbingcosts;
    float screeningcosts;
    float treatmentcosts;

    float totalcosts;
    float healthbenefitsinICU;


	//void counting_infections(); //works out who to infect
	void process_infections(gsl_rng *rng);
	void process_movements (gsl_rng *rng);
    void clear_future_events(); //clears all scheduled future movement events (otherwise run_n_sims runs into memory problems) and empty iso_queue as will

	/** moved to public for now

	void perform_events(gsl_rng *rng); //performs all scheduled events at the end of a time step
    */

}; // end of 'patient' class


