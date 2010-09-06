
/*****
HCAI_Com_IBM.cpp

Notes: just_means10a  fixes bugs with must_means10 that caused problems with the decolonization strategy (most importantly, numerical errors
 were causing pdi, progprob and pdc values to change slowly over time).

 Not sure why it's called just_means, as the code should be good for outputting everything we need shouldn't it?
 Julie? Any idea?

May  18 2010
Authors: Julie Robotham & Ben Cooper
 *
 ****/


#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <list>
#include <math.h>
#include <fstream>
#include <istream>
#include <iostream>

#include "HCAI_COM_IBM.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include </usr/local/include/gsl/gsl_randist.h>
//#include </usr/local/include/gsl/gsl_rng.h>


//#include </opt/local/include/gsl/gsl_rng.h>
//#include </opt/local/include/gsl/gsl_randist.h>

// use this on the mac include </opt/local/include/gsl/gsl_rng.h>
// use on the mac   #include </opt/local/include/gsl/gsl_randist.h>

/////////////////////CHANGED INCLUDE BIT ***************************************
//////////////////////////////////////////////////////////////////////////////////



//#include </usr/local/include/gsl/gsl_rng.h>
//#include </usr/local/include/gsl/gsl_randist.h>


//#include </usr/local/include/gsl/gsl_rng.h>


///////////////////// /////////////////////////////////////////////////////////////


using namespace std;
//extern double sensitivity;
//extern double specificity;

//code from when we were assigning a pre-set LOS on admission
//extern long unsigned int LoS;
//extern long unsigned int add_LoS;
//extern vector<unsigned long int> LoS_vector;
//extern vector<unsigned long int> add_LoS_vector;


extern bool verbose; //if true print more output

extern    int ISOCAP;

extern double effect_of_ISO;
extern double effect_of_ISO_SD;
extern double effect_of_secISO;
extern double effect_of_secISO_SD;

extern double DecolEffectonPdc_MEAN;
extern double DecolEffectonPdc_SD;
extern double DecolEffectonPdi_MEAN;
extern double DecolEffectonPdi_SD;
extern double DecolEffectonIQ_MEAN;
extern double DecolEffectonIQ_SD;
extern double DecolEffectonProgProb_MEAN;
extern double DecolEffectonProgProb_SD;

extern double Pdc;
extern double Pdi;
extern double ProgProb;
extern double IQ;
extern float prop_high_risk;

extern  double PropTreatmentSuccessful;
extern  int length_of_treatment;
//extern  float DecolEffectonPdc;
//extern  float DecolEffectonIQ;
//extern  float DecolEffectonProgProb;
extern float CC_SENSITIVITY_MEAN; //conventional culture
extern float CC_SENSITIVITY_SD; //conventional culture
extern float CC_SPECIFICITY_MEAN;
extern float CC_SPECIFICITY_SD;
extern float CA_SENSITIVITY_MEAN; //chromagar after full incubation
extern float CA_SENSITIVITY_SD; //chromagar after full incubation
extern float CA_SPECIFICITY_MEAN;
extern float CA_SPECIFICITY_SD;
extern float CA_EARLY_SENSITIVITY_MEAN; //chromagar after full incubation
extern float CA_EARLY_SENSITIVITY_SD; //chromagar after full incubation
extern float CA_EARLY_SPECIFICITY_MEAN;
extern float CA_EARLY_SPECIFICITY_SD;
extern float PCR_SENSITIVITY_MEAN;  //pcr
extern float PCR_SENSITIVITY_SD;  //pcr
extern float PCR_SPECIFICITY_MEAN;
extern float PCR_SPECIFICITY_SD;
extern    int MEANCOLDURATION;
extern    int MEANINFDURATION;
extern    int MEANHOMEDURATION;
extern    int homeduration;
extern vector<double> sens_vector;
extern vector<double> spec_vector;
extern vector<unsigned long int> age_vector;
//extern int IsoCap; //max isolation capacity


//extern vector<float> re_Prob;
extern vector<float> dis_prob_SUS;
extern vector<float> death_prob_SUS;
extern vector<float> dis_prob_INF;
extern vector<float> death_prob_INF;
extern vector<float> Pdc_vector;
extern vector<float> Pdi_vector;
extern vector<float> Progprob_vector;

extern vector<float> home_duration_SUS;


//not  using these age characteristics
extern unsigned long int age;
extern double proba;
extern double probb;
extern double probc;
extern double probd;
extern double probe;
//
extern double a; //total number of individuals in age group 0-19 in age distribution (from file)
extern double b; //total number of individuals in age group 20-39 in age distribution (from file)
extern double c; //total number of individuals in age group 40-59 age distribution (from file)
extern double d; //total number of individuals in age group 60-79 in age distribution (from file)
extern double e; //total number of individuals in age group 80+ in age distribution (from file)

extern int numsim;

extern float prop_C_on_ad;
extern float prop_C_on_ad_high_risk;
//extern  vector<screeningpolicycomponent> weekdayscreening[7]; // screening policies for each day of week
//extern  vector<screeningpolicycomponent> weeklypostadmission; // screening policies for screens every wk after a patient is admitted
//extern  vector<screeningpolicycomponent> dischargescreening;





patients::patients(gsl_rng *rng): timestepsperday(DAILYTIMESTEPS), popsize(DEFAULT_POPSIZE){ //Initiliase Patients
//IMPLEMENTATION OF THE CONSTRUCTOR FOR PATIENTS
   // int time=0; //start at time 0
    float urn; //uniform random number
    // setpolicy(rng);  //set current intervention and screening policy and parameters - set from command line, user interface or env varss
                //where command line overrides everything
  //printpolicy(); //to confirm that it was correctly set
  cumulative_StoC=0;
  cumulative_StoI=0;
  cumulative_CtoI=0;
  cumulative_Conadmission=0;
  cumulative_admissions=0;
  cumulative_readmissions=0;
  cumulative_discharges=0;
  cumulative_deaths=0;
  cumulative_appisodays=0;
  cumulative_inappisodays=0;
  cumulative_unisodays =0;
  cumulative_good_bd = 0;


if (DEFAULT_DEBUG!=16){
    //cout<<"Initialising people..."<<"\n";
}

   for(unsigned long int i=0; i<popsize; ++i){

        hos_pop[i].patientid=i;
        hos_pop[i].age=calculate_characteristics(rng);
        //cout<<hos_pop[i].age<<" chosen age for patient "<<i<<"\n";
        hos_pop[i].gender=1;
        hos_pop[i].timecolonized = 999999999;
        hos_pop[i].timeinfected = 999999999;
        hos_pop[i].timerecovered = 99999999;
        hos_pop[i].everpositiveswabthisadmission=false;
        hos_pop[i].Pdc= Pdc;//0.15;//low prev setting - calculated by trial and error to match literature scenarios//0.15 - high prev setting
        hos_pop[i].Pdi= Pdi;//0.0035;//low prev setting - calculated by trial and error to match literature scenarios//0.035;- high prev setting
	   hos_pop[i].baselinePdc=hos_pop[i].Pdc;
	   hos_pop[i].baselinePdi=hos_pop[i].Pdi;

 //       hos_pop[i].Pdc= 0.0006;//low prev setting - calculated by trial and error to match literature scenarios//0.15 - high prev setting
        //hos_pop[i].Pdi= 0.0001;//low prev setting - calculated by trial and error to match literature scenarios//0.035;- high prev setting
        hos_pop[i].ProgProb = ProgProb; //0.045;//this will need to be changed if we are assuming all C->I movements are progressions
	   hos_pop[i].baselineProgProb =hos_pop[i].ProgProb;
	   //BUT Pdi will need to have a positive value so that S->I transmission occur - in which case a C->I transmission can occur using Pdi
        hos_pop[i].IQ= IQ;//1;//get_IQ(rng);//ignore heterogeniety for now
       hos_pop[i].baselineIQ= hos_pop[i].IQ;
	   hos_pop[i].isolation=false;
        hos_pop[i].decolonization=false;
 //       hos_pop[i].decolonization_mup=false;
 //       hos_pop[i].decolonization_mup_chx=false;
        hos_pop[i].nursecohorting=false;
        hos_pop[i].previous_screen_result=NA; // i.e. no relevant previosu screen result if newly admitted
        hos_pop[i].disease_state=  SUSCEPTIBLE;
        hos_pop[i].awareness_state=UNKNOWN;
        hos_pop[i].discharge_state=NONE;
        hos_pop[i].ICU_no=99;//(99 just means not in the ICU)
        hos_pop[i].no_days_in_hos=-99; //count, while they are in ICU (admission day = 0)
        hos_pop[i].no_hos=0; //count if ever in hospital(ICU)
       // hos_pop[i].no_days_at_home=0; //number of days after discharge
        hos_pop[i].num_consec_neg_screens_following_a_pos=0; //add one for each -ve screen if prev pos (reset to 0 on a pos screen)
        hos_pop[i].specialty = GM;//at the moment this is either ICU or GM
        hos_pop[i].risk_group=LOW;
        //hos_pop[i].extra_stay_to_be_added =0; //a bit miessy having this flag as a part of the all_patients data structure

    }//end of for


    /**
    //very simply setting numbers of sus, col and inf in whole population - will change when have data
    for(unsigned long int i=0; i<10; ++i){
            hos_pop[i].disease_state=INFECTED;
    }

    for(unsigned long int i=10; i<100; ++i){
            hos_pop[i].disease_state=COLONIZED;
    }


    for(unsigned long int i=100; i<popsize; ++i){
            hos_pop[i].disease_state=SUSCEPTIBLE;
    }
    */

    for(unsigned long int i=0; i<popsize; ++i){
            hos_pop[i].disease_state=SUSCEPTIBLE;
    }//everyone starts off as being SUSCEPTIBLE, then there is a prob that this will be changed to COLONIZED in the admission routine


    //everyone starts off as being SUSCEPTIBLE, then there is a prob that this will be changed to COLONIZED in the admission routine


    // make  high risk individuals admitted more than once
    for(unsigned long int i=0; i<popsize; ++i){ //popsize WAS divided by 2.0o
        if( hos_pop[i].no_hos >=1){
   hos_pop[i].risk_group=HIGH;
  }else
  hos_pop[i].risk_group=LOW;
    }

    //choosing a random patient from list of all potential patients
    unsigned long int pick_patient;
    float rannum;

    //Choose who should start in ICU
    urn=gsl_rng_uniform(rng);
    int num_start_in_ICU = DEFAULT_ICU_SIZE;

    //struct current_ICU_patients starting_in_ICU= {NULL};//, 99, false};//order of members: pointer, bed_no, in_sideroom



    for (int k=0; k<num_start_in_ICU; ++k){//loop through those who are to start in the ICU


        rannum = gsl_rng_uniform(rng);//uniform random number: includes 0 but excludes 1
        pick_patient = int (popsize * rannum);//so will never get the final person as int rounds down??


        if (hos_pop[pick_patient].ICU_no==99){



            ICU_patients[k].ICU_patient_ptr = &hos_pop[pick_patient];//make the ICU_patient_ptr = the address of the paerson who was randomly chosen
            ICU_patients[k].bed_no=k;

            hos_pop[pick_patient].ICU_no=k;

            hos_pop[pick_patient].no_days_in_hos=0;

            hos_pop[pick_patient].no_hos=0;


        }


        else {
            //cout<<"in else, patient picked to start in ICU was already chosen"<<"\n";
            k = k-1;
        }

    }

//initial conditions for whole population
 S=popsize;//-100;
 C=10;//90
 I=0;//10
 CI=0;
 ISO=0;
 SECISO=0;
 D=0;




//initial bed day counts
 susbedday = 0;
 colbedday = 0;
 infbedday = 0;
 isobedday = 0;


//initial isolation bed day counts

    appisodays = 0;
    inappisodays = 0;
    unisodays =0;

//initial number on decol
    decolcount = 0;

//initial number in age groups
//a = 0;
//b = 0;
//c = 0;
//d = 0;
//e = 0;

//initial number of additional bed days due to infection
 addaysduetoinf=0;

//initial number of screens taken
no_ad_screens=0;
no_wkly_screens=0;
no_clin_screens=0;
no_pos_screens=0;
no_neg_screens=0;

no_pos_CC_screens=0;
no_neg_CC_screens=0;
no_pos_CA_screens=0;
no_neg_CA_screens=0;
no_pos_CA_early_screens=0;
no_neg_CA_early_screens=0;
no_PCR_screens=0;
no_ideal_screens=0;

//initial number of discharged ALIVEs and dischrage DEADs

dis_alive = 0;
dis_dead = 0;


//ICU summary variables
SinICU=DEFAULT_ICU_SIZE;
CinICU=0;
IinICU=0;
ISOinICU=0;
noICU=0;



//initial costs

beddaycosts =0;
isolationcosts=0;
decolcosts=0;
swabbingcosts = 0;
screeningcosts=0;
treatmentcosts=0;

totalcosts=0;
healthbenefitsinICU=0;


iso_queue.clear();
admis_queue.clear();


//@@@ next read in intervention parameterss - and initialize screeing parameters

} //end of initialization of 'patients' class

//IMPLEMENTATION OF DESTRUCTOR FOR PATIENTS
//patients::~patients()
//{
//}



/**printing out implemetations*/


 patients::~patients(void){         //destructor for patients
    iso_queue.clear();
    admis_queue.clear();
    future_movement_events.clear();   //need to add more here
    future_infection_events.clear();
    future_intervention_events.clear();

  }



void patients::clear_future_events(void){

 iso_queue.clear();
  admis_queue.clear();
    future_movement_events.clear();   //need to add more here
    future_infection_events.clear();
    future_intervention_events.clear();
}

 void patients::print_all_data(){

// char results_timesteps[]= "results_timesteps.txt";

 //std::ofstream resultsStream(results_timesteps);

 //resultsStream <<"Printing all data..."<<"\n";
  for (int i=0; i<DEFAULT_ICU_SIZE; ++i){
 //resultsStream <<"Age:"<<hos_pop[i].age<<"\n";
 // resultsStream <<"ICU_no:"<<hos_pop[i].ICU_no<<"\n";
//  resultsStream <<"bed_no"<<hos_pop[i].bed_no<<"\n";
  }
 }


 T_DisState patients::get_disease_state(const long int patientid){ //gets person's disease state
  return(hos_pop[patientid].disease_state);
 }



    void patients::print_all_patients_data (){
      //      cout<<"\nDisease states for all patients:  \n";
      cout<<"\nICU number for all patients:  \n";
      for(int i=0; i<DEFAULT_POPSIZE;++i){
 //cout<<"Patient "<<i<<" "
  cout    <<int(hos_pop[i].disease_state)<<" ";
 //cout    <<hos_pop[i].ICU_no<<" ";
      }
      cout<<" \n \n";
    }

 float patients::getbeddaycosts(){
  return(susbedday+colbedday+infbedday)*costofbedday;

 }

    float patients::getisolationcosts(){
        return(isobedday * dailycostofcontactprec);

 }

 float patients::getdecolcosts(){
  return(decolcount * costofdecol);

 };

 float patients::getswabbingcosts(){
  return((no_ad_screens+no_wkly_screens+no_clin_screens) * costofswab);
 };

 float patients::getscreeningcosts(){
  return((no_pos_CC_screens*costofposCC)+(no_neg_CC_screens*costofnegCC)+(no_pos_CA_screens*costofposCA)
      +(no_neg_CA_screens*costofnegCA)+(no_pos_CA_early_screens*costofposCA)+(no_neg_CA_early_screens*costofnegCA)+
      (no_PCR_screens*costofPCR));
 }

float patients::gettreatmentcosts(){
 return((cumulative_StoI + cumulative_CtoI)*costofinfetiontreatment);
}


float patients::gethealthbenefitsinICU(){
 return((susbedday+colbedday+infbedday) * (QoLinICU/365));
}

int  patients::get_col_duration(gsl_rng * rng){

         return(MEANCOLDURATION);
       //rng not used for the moment, but leave it there so we can modify to produce  an  exponentially distriubred number with mean MEANCOLDURATION
    }

    int  patients::get_inf_duration(gsl_rng * rng){
       return(MEANINFDURATION);
       //rng not used for the moment, but leave it there so we can modify to produce  an  exponentially distriubred number with mean MEANINFDURATION
    }

   int  patients::get_home_duration(gsl_rng * rng){


      /* int   prophomeduration=gsl_ran_exponential (rng, MEANHOMEDURATION);

        unsigned long int homeduration=4380-prophomeduration;

        if (homeduration>4380){

            homeduration=MEANHOMEDURATION;

            }
            else homeduration=homeduration;

            return (homeduration);

      /*}
       else
       {

         int  homeduration=1825;

        return (homeduration);*/
        return(MEANHOMEDURATION);
//       rng not used for the moment, but leave it there so we can modify to produce  an  exponentially distriubred number with mean MEANCOLDURATION
      // }

   }

/**
    int patients::get_col_duration(gsl_rng * rng){
      //returns duration of colonization samled from an exponential distriubtion in units of timesteps
     //MEANCOLDURATION must be specified in units of timesteps (can be non-interger)
     double x;
     int exponentially_distributed_duration;
     x=gsl_rng_uniform(rng);
     exponentially_distributed_duration=floor(MEANCOLDURATION*log(1/x)+0.5); //floor(x+0.5) rounds to neareset integer
     if(exponentially_distributed_duration==0) ++exponentially_distributed_duration;
     return(exponentially_distributed_duration); //above transformation produces exponentially distriubred number with mean MEANCODURATION
    }

    int patients::get_inf_duration(gsl_rng * rng){
      //returns duration of infection samled from an exponential distriubtion in units of timesteps
     //MEANINFDURATION must be specified in units of timesteps (2n be non-interger)
     double x;
     int exponentially_distributed_duration;
     x=gsl_rng_uniform(rng);
     exponentially_distributed_duration=floor(MEANINFDURATION*log(1/x)+0.5); //floor(x+0.5) rounds to neareset integer
     if(exponentially_distributed_duration==0) ++exponentially_distributed_duration;
     return(exponentially_distributed_duration); //above transformation produces exponentially distriubred number with mean MEANINFDURATION


    }


*/


             /* void patients::print_ICU_data (){//(std::vector<current_ICU_patients> ICU_patients){
  //Ben - is this code for testing & sample outputs? presumably not necessary for every simulation
    //                  to save to file




                for (int j =0; j< DEFAULT_ICU_SIZE; ++j) {//20 should be size of ICU_patients


                noICU++;




                //if (ICU_patients[j].ICU_patient_ptr->disease_state==SUSCEPTIBLE|ICU_patients[j].ICU_patient_ptr->disease_state==RECOVEREDFROMCOL|ICU_patients[j].ICU_patient_ptr->disease_state==RECOVEREDFROMINF){
                if (ICU_patients[j].ICU_patient_ptr->disease_state==SUSCEPTIBLE){
                 SinICU++;
                }

                if (ICU_patients[j].ICU_patient_ptr->disease_state==COLONIZED){
                CinICU++;
                }

//                if (ICU_patients[j].ICU_patient_ptr->disease_state==INFECTED|ICU_patients[j].ICU_patient_ptr->disease_state==INFECTEDFROMCOL){
                if (ICU_patients[j].ICU_patient_ptr->disease_state==INFECTED){
                IinICU++;
                }

                if (ICU_patients[j].ICU_patient_ptr->isolation==true){
                ISOinICU++;
                }

                cout<<" Day "<<time<<" noICU "<<noICU<<" SinICU "<<SinICU<<" CinICU "<<CinICU<<" IinICU "<<IinICU<<" ISOinICU "<<ISOinICU<<" ";
            }*/


    /* ofstream ICU_prev_file("prevalence_in_ICU.txt",  ios::app );
                ICU_prev_file<<numsim<<","<<time<<","<<SinICU<<","<<CinICU<<","<<IinICU<<",\n";

                cout<<" Day "<<time<<" SinICU "<<SinICU<<" CinICU "<<CinICU<<" IinICU "<<IinICU<<" ISOinICU "<<ISOinICU<<" ";


                //cout<<ISO<<"iso"<<"\n";

                ICU_prev_file.close();*/



                      //to print to screen
                        /////////////////////

    //    cout<<"Printing current ICU data..."<<"\n";

        /**
        for (int j =0; j< DEFAULT_ICU_SIZE; ++j) {//20 should be size of ICU_patients
            //cout<<"Age:"<<ICU_patients[j].ICU_patient_ptr->age <<"    ICU_no:"<<ICU_patients[j].ICU_patient_ptr->ICU_no<<"    bed_no:"<<ICU_patients[j].bed_no<< "\n";
            cout<<"Disease status:"<<int(ICU_patients[j].ICU_patient_ptr->disease_state) <<"    ICU_no:"<<ICU_patients[j].ICU_patient_ptr->ICU_no<<"    bed_no:"<<ICU_patients[j].bed_no;
            cout<<"  Pateint ID:"<<int(ICU_patients[j].ICU_patient_ptr->patientid)<<"\n";

            //if (ICU_patients[j].ICU_patient_ptr->interventions==ISOLATION){
            //cout<<"ID of isolated patient:"<<int(ICU_patients[j].ICU_patient_ptr->patientid)<<"\n";
            //}
        }


        */
    //}

/**MOVEMENTS*/////////////////////////////////////////////////////////////////
    void patients::process_ICU_discharges(all_patients* patient_ptr, gsl_rng * rng){




        //routine to move the person who the discharge event is happening to out of the ICU
        //(make their ICU pointer= NULL, change their ICU number and specialty code)
        //also moves them out of the IW/sideroom (if they were in it/one)

        //cout<<"doing discharges"<<"\n";


        //cout<<ISO<<"number in ISO at beginning of discharge routine"<<"\n";

        //set up movement definitions
        int moving_patient_ICU_num;
        moving_patient_ICU_num = patient_ptr->ICU_no;
        int moving_patient_ID = patient_ptr->patientid;

        struct events newevent;


        newevent.event_type=READMISSION;
        newevent.patient_ptr = patient_ptr;


        unsigned long int home_duration ;//10;//4; //currently a mean
        home_duration =get_home_duration(rng);
        duration_list[patient_ptr->patientid]=(home_duration);

       /* if (DEFAULT_DEBUG==75){

        //  cout<<home_duration<<"\t"<<"home duration"<<" \n";
//


         cout << patient_ptr->patientid <<  " \t"<< duration_list[patient_ptr->patientid] << endl;
         cout << "Vector size: " <<   " \t"<<duration_list.size() << endl;
        }*/

        //checks
//       // for (int c=1; c<32; c+=1){
//            if (DEFAULT_DEBUG==c){
//                if (moving_patient_ICU_num >DEFAULT_ICU_SIZE-1){
//                  cout<<"<moving a patient whose ICU_num exceeds ICU size"<<"\n";
//                }
//            }
//        }






        if (DEFAULT_DEBUG==3){
            cout<<int(patient_ptr->disease_state)<<"discharged patient disease state"<<"\n";
            cout<<int(patient_ptr->no_days_in_hos)<<"discharged patient LoS"<<"\n";
            cout<<int(patient_ptr->no_hos)<<"discharged patient hospitalisations"<<"\n";
            cout<<moving_patient_ICU_num<<"  discharge PATIENT ICU NUMBER"<<"\n";
            cout<<moving_patient_ID<<"   DISCHARGE PATIENT ID NUMBER " <<"\n";
        }


/**if (patient_ptr->isolation==true){
cout<<"discharged patient  id number    "<<patient_ptr->patientid<<"\n";

}
*/

        //the actual discharge
        ICU_patients[moving_patient_ICU_num].ICU_patient_ptr=NULL;
        patient_ptr->ICU_no=99;
        patient_ptr->specialty = COM;
        patient_ptr->no_days_in_hos=-99;
        patient_ptr->no_days_at_home=0; //re-set their los at home count for this episode
        patient_ptr->discharge_state = ALIVE; //discharged alive - to keep track of their QALYs post-discharge
        patient_ptr->no_hos+=1;               //count this hospitalisation
        patient_ptr->everpositiveswabthisadmission=false;
        patient_ptr->risk_group=HIGH;

       future_movement_events[time+home_duration].push_back(newevent);


    if (patient_ptr->disease_state==COLONIZED){


              --C;   //this assumes an infected patient will be replaced by a susceptible
              ++S;

            if (DEFAULT_DEBUG==72){
                cout<<int(patient_ptr->disease_state)<<" disease state at discharge"<<"\n";
            }


            ++ComC;  //Counts patients discharged colonised

        }


    //    cout<<int(patient_ptr->disease_state)<<" disease state before reset"<<"\n";

//cout<<int(patient_ptr->disease_state)<<" disease state before reset"<<"\n";
//        if (patient_ptr->disease_state==INFECTED|patient_ptr->disease_state==INFECTEDFROMCOL){
         if (patient_ptr->disease_state==INFECTED){

              //patient_ptr->disease_state=SUSCEPTIBLE;
              --I;   //this assumes an infected patient will be replaced by a susceptible
              ++S;
        }



        /** REMINDER: if there is a count associated with awareness state then this needs to be altered here!*/
        patient_ptr->awareness_state=UNKNOWN;

        //if the person discharged was in isolation, take them out
        if (patient_ptr->isolation==true){
            if (DEFAULT_DEBUG==8){
                cout<<time<<"   the person discharged from ICU was in isolation, so have been moved out " <<"\n";
            }

            patient_ptr->isolation=false;
            patient_ptr->IQ=IQ;// ignore heterogenieties for now            IQ=get_IQ(rng);//put IQ back to baseline
            --ISO;
            //cout<<patient_ptr->patientid<<"ID of patient taken out of iso due to discharge" <<"\n";
        }


        if (patient_ptr->secndisolation==true){

            if (DEFAULT_DEBUG==8){
                cout<<time<<" \t"<<" the person discharged from  ICU was in secondary isolation, so have been moved out" <<"\n";
            }

            patient_ptr->secndisolation=false;
            patient_ptr->IQ=IQ;
            --SECISO;
        }


        if (DEFAULT_DEBUG==8){
            cout<<time<<" \t"<<ISO<<" value of ISO count at end of discharge loop"<<"\n";
        }


        //if a person in the iso_queue is discharged take them off queue and stop their secondary isolation
        for (unsigned long int q=0; q<iso_queue.size(); q++){

            if (iso_queue[q]==patient_ptr->patientid){
                iso_queue.erase(iso_queue.begin()+q);

                patient_ptr->secndisolation=false;
                patient_ptr->IQ=IQ;
                --SECISO;


                if (DEFAULT_DEBUG==8){
                    cout<<*iso_queue.begin()+q<<"person taken off iso Q because they were discahrged from ICU"<<"\n";
                    cout<<"REMOVED FROM ISO_QUEUE >>>>>>>>>>>>>>>>>1"<<"\n";
                    cout<<"Patient discharged and removed from iso_queue ID number:"<<iso_queue[q]<<"\n";
                }
            }
         }


        //if discharges have created space in the IW
        //and if there are people waiting (i.e. in iso_queue)
        //move them in (but stop their sencondary isolation)
            if (ISO<policy.primary_isolation_capacity  && iso_queue.size()>0){


            if (DEFAULT_DEBUG==8){
                cout<<*iso_queue.begin()<<"person being moved from Q into ISO"<<"\n";
            }

           hos_pop[*iso_queue.begin()].isolation=true;
           ISO++;

            patient_ptr->secndisolation=false;
            --SECISO;


     //           hos_pop[*iso_queue.begin()].IQ = IQ * (1.0-effect_of_ISO);
     hos_pop[*iso_queue.begin()].IQ = IQ * (1.0-policy.primary_isolation_effectiveness);


           iso_queue.erase(iso_queue.begin());

            if (DEFAULT_DEBUG==8){
                cout<<"REMOVED FROM ISO_QUEUE >>>>>>>>>>>>>>>>>2"<<"\n";
            }
        }

        //cout<<ISO<<"number in ISO at end of discharge routine"<<"\n";

		//stop decolonization
        if (patient_ptr->decolonization==true){ //i.e. only do this if the patient is in the middle of undergoing treatment when they die
	        //cout<<"xxxx ending decol for a discharge";
			end_decolonization(patient_ptr, rng);
        }


        if(policy.dischargescreening.begin()!=policy.dischargescreening.end()) {//i.e. provided there is an discharge screening policy
                            discharge_screening(patient_ptr, rng); //perform discharge screening and schedule results of screening


                         if (DEFAULT_DEBUG==6){
                cout<<"PERFORM DISCHARGE SCREENING >>>>>>>>>>>>>>>>>"<<"\n";
            }

                        }











    } //end of process_ICU_discharges


/**

    double patients::getsensitivity(screeningtype scrtype){
     //returns senstivity of  the screening type
        switch(scrtype){

            case CC:  //conventional culture
              return(conventional_culture.sensitivity);

            break;

            case CA:  //chromagar
              return(chromagar.sensitivity);

            break;

            case CA_early:  //chromagar early result
              return(chromagar_early.sensitivity);

            break;

            case PCR:
              return(pcr.sensitivity);

            break;

            default:
                    cout<<"\n Error: Undefined screening type used\n";
                    exit(0);

        }
    }

     double patients::getspecificity(screeningtype scrtype){
     //returns specificity of  the screening type
        switch(scrtype){

            case CC:  //conventional culture
              return(conventional_culture.specificity);

            break;

            case CA:  //chromagar
              return(chromagar.specificity);

            break;

            case CA_early:  //chromagar early result
              return(chromagar_early.specificity);

            break;

            case PCR:
              return(pcr.specificity);

            break;

            default:
                    cout<<"\n Error: Undefined screening type used\n";
                    exit(0);

        }

     }

     double patients::gettat(screeningtype scrtype){
    //returns turnaround time of  the screening type
        switch(scrtype){

            case CC:  //conventional culture
              return(conventional_culture.tat);

            break;

            case CA:  //chromagar
              return(chromagar.tat);

            break;

            case CA_early:  //chromagar early result
              return(chromagar_early.tat);

            break;

            case PCR:
              return(pcr.tat);

            break;

            default:
                    cout<<"\n Error: Undefined screening type used\n";
                    exit(0);

        }
     }

*/

    void patients::process_ICU_deaths(all_patients* patient_ptr, gsl_rng * rng){
        //moved deaths out of ICU
        //(make their ICU pointer= NULL, change their ICU number and specialty code)
        //also moves them out of the IW/sideroom (if they were in it/one)






        int dying_patient_ICU_num;
        dying_patient_ICU_num = patient_ptr->ICU_no;

      /*  for (int c=1; c<32; c+=1){
            if (DEFAULT_DEBUG==c){
                if (dying_patient_ICU_num >19){
                cout<<"<STOP!!!!!!!!!!!THIS ISNT ACTUALLY WORKING!"<<"\n";
                }
            }
        }*/

       // int dying_patient_ID = patient_ptr->patientid;


                //cout<<dying_patient_ICU_num<<"  DYING PATIENT ICU NUMBER"<<"\n";
                //cout<<dying_patient_ID<<"   DYING PATIENT ID NUMBER " <<"\n";


        ICU_patients[dying_patient_ICU_num].ICU_patient_ptr=NULL;

        patient_ptr->ICU_no=99;
        patient_ptr->specialty = NA;
        patient_ptr->no_days_in_hos=-99;
        patient_ptr->no_days_at_home=0;//reset their los count for thsi hospital episode
        patient_ptr->discharge_state = DEAD;


        if (patient_ptr->disease_state==COLONIZED){

            if (DEFAULT_DEBUG==1){
                cout<<int(patient_ptr->disease_state)<<" disease state before reset"<<"\n";
            }

            patient_ptr->disease_state=SUSCEPTIBLE;
            --C;
            ++S;
        }

     //   if (patient_ptr->disease_state==INFECTED|patient_ptr->disease_state==INFECTEDFROMCOL){
        if (patient_ptr->disease_state==INFECTED){
            if (DEFAULT_DEBUG==1){
            //cout<<int(patient_ptr->disease_state)<<" disease state before reset"<<"\n";
            }
              patient_ptr->disease_state=SUSCEPTIBLE;
              --I;    //default assumption is that infected is replaced by as suceptible
              ++S;
        }



        //if the person who died was in isolation, take them out
        if (patient_ptr->isolation==true){

            if (DEFAULT_DEBUG==8){
                cout<<"the person who died in ICU was in isolation, so have been moved out" <<"\n";
            }

            patient_ptr->isolation=false;
            patient_ptr->IQ=IQ;
            --ISO;

            //cout<<patient_ptr->patientid<<"ID of patient taken out of iso due to death" <<"\n";
        }


       if (patient_ptr->secndisolation==true){

            if (DEFAULT_DEBUG==8){
                cout<<"the person who died in ICU was in secondary isolation, so have been moved out" <<"\n";
            }

            patient_ptr->secndisolation=false;
            patient_ptr->IQ=IQ;
            --SECISO;
        }


        if (DEFAULT_DEBUG==8){
            cout<<ISO<<"value of ISO count at end of death loop"<<"\n";
            cout<<SECISO<<"value of SECndISO count at end of death loop"<<"\n";
        }

        //if a person in the iso_queue dies take them off queue and stop secondary isolation for them
        for (unsigned long int q=0; q<iso_queue.size(); q++){
            if (iso_queue[q]==patient_ptr->patientid){

                iso_queue.erase(iso_queue.begin()+q);
                patient_ptr->secndisolation=false;
                patient_ptr->IQ=IQ;
                --SECISO;

                if (DEFAULT_DEBUG==8){
                    cout<<*iso_queue.begin()+q<<"person taken off iso Q because they died"<<"\n";

                    //cout<<"REMOVED FROM ISO_QUEUE DUE TO DEATH>>>>>>>>>>>>>>>>>1"<<"\n";
                    //cout<<"Patient died and removed from iso_queue ID number:"<<iso_queue[q]<<"\n";

                }

            }
         }


        //if deaths have created space in the IW
        //and if there are people waiting (i.e. in iso_queue)
        //move them in and stop their secondary isolation
        if (ISO<policy.primary_isolation_capacity   && iso_queue.size()>0){

                if (DEFAULT_DEBUG==8){
                    cout<<*iso_queue.begin()<<"person being moved from Q into ISO"<<"\n";
                }

               hos_pop[*iso_queue.begin()].isolation=true;
               ISO++;
               patient_ptr->secndisolation=false;
               --SECISO;

        //               hos_pop[*iso_queue.begin()].IQ = IQ * (1-effect_of_ISO);
               hos_pop[*iso_queue.begin()].IQ = IQ * (1-policy.primary_isolation_effectiveness);


               iso_queue.erase(iso_queue.begin());

                if (DEFAULT_DEBUG==8){
                    cout<<"REMOVED FROM ISO_QUEUE >>>>>>>>>>>>>>>>>2"<<"\n";
                }
        }

		//stop decolonization
        if (patient_ptr->decolonization==true){ //i.e. only do this if the patient is in the middle of undergoing treatment when they die
		  //cout<<"xxxx ending decol for a death";
			end_decolonization(patient_ptr, rng);
        }


    }//end of process_ICU_deaths




    void patients::process_ICU_admissions(gsl_rng * rng){

        //routine to pick a random patient from from whole pop,
        //and move them into an empty bed in the ICU

        //cout<<"in process ICU admissions"<<"\n";

        unsigned long int pick_patient;
        float rannum;
        rannum = gsl_rng_uniform(rng);//uniform random number: includes 0 but excludes 1
        pick_patient = int (popsize * rannum);//so will never get the final person as int rounds down???????

        //cout<<hos_pop[pick_patient].ICU_no<<" ICU number of paicked patient, if less than 19 need to pick someone else" <<"\n";

 //    print_all_patients_data();
 // int xxx; std::cin>>xxx;
 //int zz=0;

//patients not currently in the ICU should all have an ICU number of 99, while this in the ICU have an ICU number that equals their
// bed number ... so to make sure we don't pick someone not already in the icu, keep trying while ICU_no is not equal to 99


    while(hos_pop[pick_patient].ICU_no!=99){ //i.e. keep trying to pick a patient until we find one not in ICU (so ICU_no==99)
        rannum = gsl_rng_uniform(rng);//uniform random number: includes 0 but excludes 1
        pick_patient = int (popsize * rannum);

 }

        int admitC = 0;

        if (DEFAULT_DEBUG==3){

            cout<<hos_pop[pick_patient].no_hos<<" number of hosp of chosen admittee" <<"\n";

             cout<<int(hos_pop[pick_patient].disease_state)<<" disease state of chosen admittee"<<"\n";

        }





if(hos_pop[pick_patient].decolonization) cout<<"Error! xxxxxxx. Admitted patient already undergoing decolonization";

 //-------------------set disease state of admitted person-----------------------------------------------------

     float random;
     random = gsl_rng_uniform(rng);
     bool makecolonized=false;


      //Risk status is set at discharge currently




    if (hos_pop[pick_patient].no_hos>=1){

           hos_pop[pick_patient].risk_group=HIGH;
           }
           else{
           hos_pop[pick_patient].risk_group=LOW;
           }


if (hos_pop[pick_patient].risk_group==HIGH){
	if (random<prop_C_on_ad_high_risk){  //then make colonized
    makecolonized=true;
    } else {
    //  do nothing
    }
   } else { //not high risk

     if (random<prop_C_on_ad){  //draw random number to see if they are Colonized on admission
      makecolonized=true;
     }
   }
   if(makecolonized){  //then default assumption that new admission is suceptible is wrong, so we increment C and decrement S
                 //    if (hos_pop[pick_patient].disease_state==SUSCEPTIBLE|hos_pop[pick_patient].disease_state==RECOVEREDFROMINF|hos_pop[pick_patient].disease_state==RECOVEREDFROMCOL){
                    if (hos_pop[pick_patient].disease_state==SUSCEPTIBLE){
                        S--;
                       C++;
                        hos_pop[pick_patient].disease_state=COLONIZED;




                    }
   }



 //hos_pop[pick_patient].disease_state=COLONIZED;

            if (DEFAULT_DEBUG==3){
                    cout<<int(hos_pop[pick_patient].disease_state)<<" disease state of chosen admittee"<<"\n";
            }

            if (DEFAULT_DEBUG==3){
            cout<<" Number admitted colonised"<<" \t"<<admitC<<" \n";

            }




   if (DEFAULT_DEBUG==3){
                    cout<<int(hos_pop[pick_patient].disease_state)<<" disease state of chosen admittee"<<"\n";
            }

            if (DEFAULT_DEBUG==3){
            cout<<" Number admitted colonised"<<" \t"<<admitC<<" \n";

            }
//______________________________CALLING SCREENING ON ADMISSION_____________________________//


                /**screening of high risk patients on admission*/
                               // if (time>START_TIME_OF_INT){
                if (time>policy.timesteptoimplement){ //i.e. if policy has been implented by current time step

                    if(policy.screeninghighriskonly) {

                       // cout<<"should only be doing this if we are screening high risk only"<<"\n";

                        if(policy.admissionscreening.begin()!=policy.admissionscreening.end()) {
                            if (hos_pop[pick_patient].risk_group==HIGH){
                            admission_screening(rng, pick_patient); //perform admission screening and schedule results of screening
                            }
                        }
                    }
                    else{


                            /**screening of patients on admission*/

                      //  cout<<"should only be doing this if we ARENT screening high risk only"<<"\n";

                        if(policy.admissionscreening.begin()!=policy.admissionscreening.end()) {//i.e. provided there is an admission screening policy
                            admission_screening(rng, pick_patient); //perform admission screening and schedule results of screening
                        }
                    }




                }

//______________________________CALLING INTERVENTIONS ON ADMISSION_____________________________//

                /**interventions on admission pre-emptive isolation and blanket decol (of all and just high risk)*/





                if (time>policy.timesteptoimplement){ //i.e. if policy has been implented by current time ste
                      implement_control_measures_for_all(rng, pick_patient);


                if (hos_pop[pick_patient].risk_group==HIGH){
                    implement_control_measures_for_high_risk(rng, pick_patient);
                }

}


                if (DEFAULT_DEBUG==1){
                    //cout<<hos_pop[pick_patient].ICU_no<<"   ICU_no before admission" <<&hos_pop[pick_patient]<<"their address"<<"\n";
                }



//choosing empty bed that admitted patient is moved in to //

                for (int m=0; m<DEFAULT_ICU_SIZE; m++){


                    if (ICU_patients[m].ICU_patient_ptr==NULL){/**what about if there is more than- 1 empty bed?*/

                            if (DEFAULT_DEBUG==1){
                                cout<<ICU_patients[m].bed_no<<"bed_no";
                            }

                            ICU_patients[m].ICU_patient_ptr=&hos_pop[pick_patient];


                            if (DEFAULT_DEBUG==1){
                                cout<<ICU_patients[m].ICU_patient_ptr->patientid<<"  PATIENT ID OF ADMITTED PERSON" <<"\n";
                            }

                            //cout<<ICU_patients[m].ICU_patient_ptr->patientid<<"  PATIENT ID OF ADMITTED PERSON" <<"\n";

                            ICU_patients[m].ICU_patient_ptr->ICU_no=m;

                            ICU_patients[m].ICU_patient_ptr->specialty=ICU;

                            ICU_patients[m].ICU_patient_ptr->no_days_in_hos=0;


                            if (DEFAULT_DEBUG==1){
                                cout<<ICU_patients[m].ICU_patient_ptr->ICU_no<<"  ICU_no after admission"<<"\n";
                                cout<<int(ICU_patients[m].ICU_patient_ptr->disease_state)<<"  disease state of admittee"<<"\n";
                            }


                    }//if bed is empty
                }//for

    //cout<<"Day "<<time<<" admitC "<<admitC<<" ";
//temporary debugging code xxxxxx may 2010
/*
		hos_pop[pick_patient].Pdc= Pdc;//0.15;//low prev setting - calculated by trial and error to match literature scenarios//0.15 - high prev setting
        hos_pop[pick_patient].Pdi= Pdi;//0.0035;//low prev setting - calculated by trial and error to match literature scenarios//0.035;- high prev setting

		//       hos_pop[i].Pdc= 0.0006;//low prev setting - calculated by trial and error to match literature scenarios//0.15 - high prev setting
        //hos_pop[i].Pdi= 0.0001;//low prev setting - calculated by trial and error to match literature scenarios//0.035;- high prev setting
        hos_pop[pick_patient].ProgProb = ProgProb; //0.045;//this will need to be changed if we are assuming all C->I movements are progressions
*/

    }// process_ICU_admissions




/**INTERVENTIONS*/////////////////////////////////////////////////////////////////

    void patients::discharge_screening(all_patients* patient_ptr, gsl_rng *rng){

    //performs admission screening
    //for each admitted patient:
    //   - chooses if they are screened
    // schedules an act on screen event which occurs at time TAT in the future
    // Note that future events include changing their awareness_state
    // note also that screening could possibly include more than one type of screening
    // or different types for different patients
if (time>policy.timesteptoimplement){

        struct events newevent1;
        struct events newevent2;
        struct events newevent3;
        newevent1.event_type=SCHEDULED_SCREEN; //i.e. if we are scheduleing a future screening
        newevent2.event_type=ACT_ON_POS_SCREEN;
        newevent3.event_type=ACT_ON_NEG_SCREEN;
        newevent1.patient_ptr = &hos_pop[patient_ptr->patientid];
        newevent2.patient_ptr = &hos_pop[patient_ptr->patientid];
        newevent3.patient_ptr = &hos_pop[patient_ptr->patientid];

  //      int rescreen_interval=7*timestepsperday; //specified in case  weekly rescreening is implemented

        float rannum, rannum1;
        bool screenthispatient; //patients are screened weekly with prob policy1.proportionscreenedweekly

        rannum1 = gsl_rng_uniform(rng);

        int tat;           //turnaround time
        float sens, spec; //sensitivity and secificity


        rannum = gsl_rng_uniform(rng);
        if(rannum<policy.proportionscreenedondischarge) screenthispatient=true; else screenthispatient=false;//check to see if we select this patient
        //set default values to those for conventional cultures
        //sens=conventional_culture.sensitivity;
        //spec=conventional_culture.specificity;
        //tat=conventional_culture.tat;

        //loop through vector of amdission screening swabs (since possibly more than one type of admissions swab
        for(vector<screeningpolicycomponent>::const_iterator iter=policy.dischargescreening.begin(); iter!=policy.dischargescreening.end();++iter){

          //if (rannum1<policy.proportionscreenedonadmission){ //to account for fact that proportion screened may be <1 if compliance not perfect
                if(iter->technology==CC) {sens=conventional_culture.sensitivity;spec=conventional_culture.specificity;tat=conventional_culture.tat;}
                if(iter->technology==CA) {sens=chromagar.sensitivity;spec=chromagar.specificity;tat=chromagar.tat;}
                if(iter->technology==CA_early) {sens=chromagar_early.sensitivity;spec=chromagar_early.specificity;tat=chromagar_early.tat;}
                if(iter->technology==PCR) {sens=pcr.sensitivity;spec=pcr.specificity;tat=pcr.tat;}
                if(iter->technology==I) {sens=ideal.sensitivity;spec=ideal.specificity;tat=ideal.tat;}


                if (DEFAULT_DEBUG==6){
                    cout<<int (hos_pop[patient_ptr->patientid].disease_state)<<"\t"<<"PATIENT ID OF DISCHARGE SCREENEE"<<"\n";
                }


                if (hos_pop[patient_ptr->patientid].disease_state==COLONIZED||hos_pop[patient_ptr->patientid].disease_state==INFECTED){

                    float rannum2;
                    rannum2 = gsl_rng_uniform(rng);
                    //cout<<sensitivity<<"sensitivity when using it"<<"\n";


                    if (rannum2<sens){ //i.e. true case is detected

                        if (DEFAULT_DEBUG==6){
                            cout<<"SCREENED ON DISCHARGE AND FOUND COL (TRUE POS)"<< "  person="<<hos_pop[patient_ptr->patientid].patientid<<"\n";
                        }
//  note that awareness state shouldn't change straight away......we should only act on it after tat
//- so need to shedule an act on a positive screen which should change the awareness state
//
                        ++no_pos_screens;
                        if(iter->technology==CC){
                            ++no_pos_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_pos_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_pos_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }


//cout<<" queued an act on positive screen (rightly) for patient "<<hos_pop[pick_patient].patientid<<"\n";
//cout<<" at time " << time+tat <<"\n";

                        future_intervention_events[time+tat].push_back(newevent2); //schedule act on positive screenscreen


                    } //if (rannum2<sensitivity)


                    else {//if colonized person is screened inaccurately (i.e. false negative screen)

                        if (DEFAULT_DEBUG==6){
                            cout<<"SCREENED AND FOUND NEG (FALSE NEG)at discharge"<<"  person="<<hos_pop[patient_ptr->patientid].patientid<<"\n";
                        }

                        ++no_neg_screens;

                        if(iter->technology==CC){
                            ++no_neg_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_neg_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_neg_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }

//cout<<" queued an act on negative screen (wrongly) for patient "<<hos_pop[pick_patient].patientid<<"\n";
//cout<<" at time " << time+tat <<"\n";



                        future_intervention_events[time+tat].push_back(newevent3); //schedule act on negative screen
                    }
                } else if(hos_pop[patient_ptr->patientid].disease_state==SUSCEPTIBLE) { //if patient is not colonized or infected they must be susceptible


                    float rannum5;
                    rannum5 = gsl_rng_uniform(rng);

                    if (rannum5>spec){
                        if (DEFAULT_DEBUG==6){
                            cout<<"FALSE POSITIVE SWAB FROM SCREENING ON DISCHARGE"<<"  person="<<hos_pop[patient_ptr->patientid].patientid<<"\n";
                        }

                        ++no_pos_screens;


                        if(iter->technology==CC){
                            ++no_pos_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_pos_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_pos_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                        ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }

//cout<<" queued an act on positive screen (wrongly) for patient "<<hos_pop[pick_patient].patientid<<"\n";
//cout<<" at time " << time+tat <<"\n";

                        future_intervention_events[time+tat].push_back(newevent2); //schedule act on positive screenscreen

                    } else {
                        if (DEFAULT_DEBUG==6){
                            cout<<"TRUE NEGATIVE SWAB"<<"  person="<<hos_pop[patient_ptr->patientid].patientid<<"\n";
                        }

                        ++no_neg_screens;


                        if(iter->technology==CC){
                            ++no_neg_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_neg_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_neg_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }

//cout<<" queued an act on negative screen (rightly) for patient "<<hos_pop[pick_patient].patientid<<"\n";
//cout<<" at time " << time+tat <<"\n";

                        future_intervention_events[time+tat].push_back(newevent3); //schedule act on negative screen
                    }
                } else {  //if patient is not colonized or infected or susceptible  232 `

                    if (DEFAULT_DEBUG==6){
                    cout<<"\n disease_state not SUSCEPTIBLE, INFECTED or COLONIZED" ;

                    }
                    exit(0);
                }
          }
        }
}
    void patients::admission_screening(gsl_rng *rng, int pick_patient){

    //performs admission screening
    //for each admitted patient:
    //   - chooses if they are screened
    // schedules an act on screen event which occurs at time TAT in the future
    // Note that future events include changing their awareness_state
    // note also that screening could possibly include more than one type of screening
    // or different types for different patients
if (time>policy.timesteptoimplement){

        struct events newevent1;
        struct events newevent2;
        struct events newevent3;
        newevent1.event_type=SCHEDULED_SCREEN; //i.e. if we are scheduleing a future screening
        newevent2.event_type=ACT_ON_POS_SCREEN;
        newevent3.event_type=ACT_ON_NEG_SCREEN;
        newevent1.patient_ptr = &hos_pop[pick_patient];
        newevent2.patient_ptr = &hos_pop[pick_patient];
        newevent3.patient_ptr = &hos_pop[pick_patient];

        int rescreen_interval=7*timestepsperday; //specified in case  weekly rescreening is implemented

        float rannum, rannum1;
        bool screenthispatient; //patients are screened weekly with prob policy1.proportionscreenedweekly

        rannum1 = gsl_rng_uniform(rng);

        int tat;           //turnaround time
        float sens, spec; //sensitivity and secificity


        rannum = gsl_rng_uniform(rng);
        if(rannum<policy.proportionscreenedonadmission) screenthispatient=true; else screenthispatient=false;//check to see if we select this patient
        //set default values to those for conventional cultures
        //sens=conventional_culture.sensitivity;
        //spec=conventional_culture.specificity;
        //tat=conventional_culture.tat;


        //SCHEDULE THEIR WEEKLY RE-SCREEN if we are doing this
        if(isthereweeklypostadmissionscreening() && screenthispatient){
            future_intervention_events[time+rescreen_interval].push_back(newevent1);
        }

        //loop through vector of amdission screening swabs (since possibly more than one type of admissions swab
        for(vector<screeningpolicycomponent>::const_iterator iter=policy.admissionscreening.begin(); iter!=policy.admissionscreening.end();++iter){

          //if (rannum1<policy.proportionscreenedonadmission){ //to account for fact that proportion screened may be <1 if compliance not perfect
                if(iter->technology==CC) {sens=conventional_culture.sensitivity;spec=conventional_culture.specificity;tat=conventional_culture.tat;}
                if(iter->technology==CA) {sens=chromagar.sensitivity;spec=chromagar.specificity;tat=chromagar.tat;}
                if(iter->technology==CA_early) {sens=chromagar_early.sensitivity;spec=chromagar_early.specificity;tat=chromagar_early.tat;}
                if(iter->technology==PCR) {sens=pcr.sensitivity;spec=pcr.specificity;tat=pcr.tat;}
                if(iter->technology==I) {sens=ideal.sensitivity;spec=ideal.specificity;tat=ideal.tat;}


                if (DEFAULT_DEBUG==2){
                    cout<<int (hos_pop[pick_patient].disease_state)<<"\n";
                }


                if (hos_pop[pick_patient].disease_state==COLONIZED||hos_pop[pick_patient].disease_state==INFECTED){

                    float rannum2;
                    rannum2 = gsl_rng_uniform(rng);
                    //cout<<sensitivity<<"sensitivity when using it"<<"\n";


                    if (rannum2<sens){ //i.e. true case is detected

                        if (DEFAULT_DEBUG==2){
                            cout<<"SCREENED AND FOUND COL (TRUE POS)"<< "  person="<<hos_pop[pick_patient].patientid<<"\n";
                        }
//  note that awareness state shouldn't change straight away......we should only act on it after tat
//- so need to shedule an act on a positive screen which should change the awareness state
//
                        ++no_pos_screens;
                        if(iter->technology==CC){
                            ++no_pos_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_pos_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_pos_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }


//cout<<" queued an act on positive screen (rightly) for patient "<<hos_pop[pick_patient].patientid<<"\n";
//cout<<" at time " << time+tat <<"\n";

                        future_intervention_events[time+tat].push_back(newevent2); //schedule act on positive screenscreen


                    } //if (rannum2<sensitivity)


                    else {//if colonized person is screened inaccurately (i.e. false negative screen)

                        if (DEFAULT_DEBUG==2){
                            cout<<"SCREENED AND FOUND NEG (FALSE NEG)"<<"  person="<<hos_pop[pick_patient].patientid<<"\n";
                        }

                        ++no_neg_screens;

                        if(iter->technology==CC){
                            ++no_neg_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_neg_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_neg_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }

//cout<<" queued an act on negative screen (wrongly) for patient "<<hos_pop[pick_patient].patientid<<"\n";
//cout<<" at time " << time+tat <<"\n";



                        future_intervention_events[time+tat].push_back(newevent3); //schedule act on negative screen
                    }
                } else if(hos_pop[pick_patient].disease_state==SUSCEPTIBLE) { //if patient is not colonized or infected they must be susceptible


                    float rannum5;
                    rannum5 = gsl_rng_uniform(rng);

                    if (rannum5>spec){
                        if (DEFAULT_DEBUG==2){
                            cout<<"FALSE POSITIVE SWAB"<<"  person="<<hos_pop[pick_patient].patientid<<"\n";
                        }

                        ++no_pos_screens;


                        if(iter->technology==CC){
                            ++no_pos_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_pos_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_pos_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                        ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }

//cout<<" queued an act on positive screen (wrongly) for patient "<<hos_pop[pick_patient].patientid<<"\n";
//cout<<" at time " << time+tat <<"\n";

                        future_intervention_events[time+tat].push_back(newevent2); //schedule act on positive screenscreen

                    } else {
                        if (DEFAULT_DEBUG==2){
                            cout<<"TRUE NEGATIVE SWAB"<<"  person="<<hos_pop[pick_patient].patientid<<"\n";
                        }

                        ++no_neg_screens;


                        if(iter->technology==CC){
                            ++no_neg_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_neg_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_neg_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }

//cout<<" queued an act on negative screen (rightly) for patient "<<hos_pop[pick_patient].patientid<<"\n";
//cout<<" at time " << time+tat <<"\n";

                        future_intervention_events[time+tat].push_back(newevent3); //schedule act on negative screen
                    }
                } else {  //if patient is not colonized or infected or susceptible  232 `
                    cout<<"\n disease_state not SUSCEPTIBLE, INFECTED or COLONIZED" ;
                    exit(0);
                }
//move to act on screen bit @@@@@@@
//@@@ but also change this so that if they had a previosu positive screen on same stay, it takes three negative screens before
//@@@ their status changes - so need to keep track of  how many negative screens previously positive people have had
//                       hos_pop[pick_patient].awareness_state = BELIEVED_S_ACTUALLY_C;

//@@@@@@ change this to schedule an ACT_ON_NEG_SCREENevent (and change ACT_ON_SCREEN to ACT_ON_POS_SCREEN
                            //if patient was in isoaltion (from pre-emptive iso then move them out based on screening results).
//                            if (hos_pop[pick_patient].interventions==ISOLATION){
//                                hos_pop[pick_patient].interventions=NONE;
//                                hos_pop[pick_patient].IQ= 1; //for baseline IQ from distribution: get_IQ(rng);
//                                --ISO;
//                            }

//@@@@@@ change stuff below to schedule an ACT_ON_NEG_SCREENevent which happen only after TAT
                            //if screens have created space in the IW
                            //and if there are people waiting (i.e. in iso_queue)
                            //move them in
//                            if (ISO<IsoCap && iso_queue.size()>0){


//                               if (DEFAULT_DEBUG==8){
//                                    cout<<*iso_queue.begin()<<"person moved out of Q into ISO because rescreens have made space "<<"\n";
 //                              }

 //                               hos_pop[*iso_queue.begin()].interventions=ISOLATION;
 //                               ISO++;
 //                               hos_pop[*iso_queue.begin()].IQ= 1 * (1-effect_of_ISO);

 //                               iso_queue.erase(iso_queue.begin());

 //                               if (DEFAULT_DEBUG==8){
 //                                   cout<<"REMOVED FROM ISO_QUEUE >>>>>>>>>>>>>>>>>7"<<"\n";
  //                              }
//
//
//                            }
//
//
//@@@@@@ change stuff below to schedule an ACT_ON_NEG_SCREENevent which happen only after TAT
//                            //if a person in the iso_queue is screened and found -ve take them off queue
//                            for (unsigned long int q=0; q<iso_queue.size(); q++){
//
//
//
//                                if (iso_queue[q]==hos_pop[pick_patient].patientid){
//
//                                    if (DEFAULT_DEBUG==8){
//                                        cout<<hos_pop[pick_patient].patientid<<"patient ID of person rescreened and found -ve compared to "<<*iso_queue.begin()+q<<"ID number of person removed from Q"<<"\n";
//                                    }
//
//                                   iso_queue.erase(iso_queue.begin()+q); //dont know if this isactually erasing the right person!!!!!
//
//                                   if (DEFAULT_DEBUG==8){
//                                        cout<<"REMOVED FROM ISO_QUEUE >>>>>>>>>>>>>>>>>8"<<"\n";
//                                   }
//
//                                }
//                             }//end of for
//
//                    }//end of else i.e. when they are screened inaccurately
//
//        }//end of if - if they are colonized

//            //screening of INFECTEDS
 //           if (hos_pop[pick_patient].disease_state==INFECTED|hos_pop[pick_patient].disease_state==INFECTEDFROMCOL){
//           if (hos_pop[pick_patient].disease_state==INFECTED){
//
//                                        //cout<<"got to here inf"<<"\n";
//
//                        float rannum3;
//                        rannum3 = gsl_rng_uniform(rng);
//                        if (rannum3<sensitivity){
//
//                            if (DEFAULT_DEBUG==2){
//                            cout<<"SCREENED AND FOUND INF"<<"  person="<<hos_pop[pick_patient].patientid<<"\n";
//                            }
//
//@@@@@@ change stuff below to schedule an ACT_ON_POS_SCREENevent which happen only after TAT
//
//                            hos_pop[pick_patient].awareness_state = KNOWN_I;
//
//
//                            //SCHEDULE THEIR RE-SCREEN
//                            @@@@ change this so that rather than having a rescreen after rescreen interval
//                            rescreen occurs at fixed time intervals
//
//                            future_intervention_events[time+rescreen_interval].push_back(newevent1);
//                            future_intervention_events[time+TAT].push_back(newevent2);
//
//                            if (DEFAULT_DEBUG==2){
//                               cout<<"Act on screen scheduled for person="<<hos_pop[pick_patient].patientid<<"\n";
//                            }
//
//                            /** CALLING INTERVENTIONS FROM SCREENING
//
//                            //PUT THEM IN ISOLATION
//                            isolation(pick_patient, rng);
//
//
//                            if (time>100){
//                            decolonization(rng, pick_patient);
//                            }
//                            */
//
//                        }
//
//                        else{//if a false negative screen
//@@@@@@ change stuff below to schedule an ACT_ON_NEG_SCREENevent which happen only after TAT
//                            if (DEFAULT_DEBUG==2){
//                            cout<<"Believed S actually I"<<"  person="<<hos_pop[pick_patient].patientid<<"\n";
//                            }
//
//                            hos_pop[pick_patient].awareness_state = BELIEVED_S_ACTUALLY_I;
//
//
//                            //if patient was in isoaltion (from pre-emptive iso then move them out based on screening results).
//                            if (hos_pop[pick_patient].interventions==ISOLATION){
//                                hos_pop[pick_patient].interventions=NONE;
//                                hos_pop[pick_patient].IQ= 1; //for baseline IQ from distribution: get_IQ(rng);
//                                --ISO;
//                            }
//
//
//                            //if screens have created space in the IW
//                            //and if there are people waiting (i.e. in iso_queue)
//                            //move them in
//                            if (ISO<IsoCap && iso_queue.size()>0){
//
//
//                               if (DEFAULT_DEBUG==8){
//                                    cout<<*iso_queue.begin()<<"person moved out of Q into ISO because rescreens have made space "<<"\n";
//                               }
//
//                                hos_pop[*iso_queue.begin()].interventions=ISOLATION;
//                                ISO++;
//                                hos_pop[*iso_queue.begin()].IQ= 1 * (1-effect_of_ISO);
//
//                                iso_queue.erase(iso_queue.begin());
//
//                                if (DEFAULT_DEBUG==8){
//                                    cout<<"REMOVED FROM ISO_QUEUE >>>>>>>>>>>>>>>>>7"<<"\n";
//                                }
//
//                            }
//
//
//                            //if a person in the iso_queue is screened and found -ve take them off queue
//                            for (unsigned long int q=0; q<iso_queue.size(); q++){
//
//
//
//                                if (iso_queue[q]==hos_pop[pick_patient].patientid){
//
//                                    if (DEFAULT_DEBUG==8){
//                                        cout<<hos_pop[pick_patient].patientid<<"patient ID of person rescreened and found -ve compared to "<<*iso_queue.begin()+q<<"ID number of person removed from Q"<<"\n";
//                                    }
//
//                                   iso_queue.erase(iso_queue.begin()+q); //dont know if this isactually erasing the right person!!!!!
//
//                                   if (DEFAULT_DEBUG==8){
//                                        cout<<"REMOVED FROM ISO_QUEUE >>>>>>>>>>>>>>>>>8"<<"\n";
//                                   }
//
//                                }
//                             }
//
//
//                        }//end of else - when infected person was screened inaccurately
//            }//end of if  - if they are infected

            ////screening of SUSCEPTIBLES

            //if (hos_pop[pick_patient].disease_state==SUSCEPTIBLE|hos_pop[pick_patient].disease_state==RECOVEREDFROMINF|hos_pop[pick_patient].disease_state==RECOVEREDFROMCOL){
//            if (hos_pop[pick_patient].disease_state==SUSCEPTIBLE){
//
//                float rannum5;
//                //cout<<"got to here"<<"\n";
//                rannum5 = gsl_rng_uniform(rng);
//                //cout<<specificity<<"\n";
//
//                if (rannum5<specificity){
//                //cout<<specificity<<"specificity"<<"\n";
//
//                    if (DEFAULT_DEBUG==2){
//                        cout<<"SCREENED AND FOUND SUS"<<"  person="<<hos_pop[pick_patient].patientid<<"\n";
//                    }
//
//                    //SCHEDULE THEIR RE-SCREEN
//                    future_intervention_events[time+rescreen_interval].push_back(newevent1);
//                    hos_pop[pick_patient].awareness_state = KNOWN_S;
//
//
//                    //if patient was in isoaltion (from pre-emptive iso then move them out based on screening results).
//@@@ this should only be done after TAT - so change this
//                    if (hos_pop[pick_patient].interventions==ISOLATION){ //so currently 1 -ve screen is enough to unisolate
//                        hos_pop[pick_patient].interventions=NONE;
//                        hos_pop[pick_patient].IQ= 1; //for baseline IQ from distribution: get_IQ(rng);
//                        --ISO;
//                    }
//
//@@ again...this can only happen after TAT...since nothing happens immediately when screen are taken
//                    //if screens have created space in the IW
//                    //and if there are people waiting (i.e. in iso_queue)
//                    //move them in
//                    if (ISO<IsoCap && iso_queue.size()>0){
//
//
//                        if (DEFAULT_DEBUG==8){
//                            cout<<*iso_queue.begin()<<"person moved out of Q into ISO because rescreens have made space "<<"\n";
//                        }
//
//                        hos_pop[*iso_queue.begin()].interventions=ISOLATION;
//                        ISO++;
//                        hos_pop[*iso_queue.begin()].IQ= 1 * (1-effect_of_ISO);
//
//                        iso_queue.erase(iso_queue.begin());
//
//                        if (DEFAULT_DEBUG==8){
//                            cout<<"REMOVED FROM ISO_QUEUE >>>>>>>>>>>>>>>>>7"<<"\n";
//                        }
//
//                    }
//
//                    //if a person in the iso_queue is screened and found -ve take them off queue
//                            for (unsigned long int q=0; q<iso_queue.size(); q++){
//
//
//                                if (iso_queue[q]==hos_pop[pick_patient].patientid){
//
//                                    if (DEFAULT_DEBUG==8){
//                                        cout<<hos_pop[pick_patient].patientid<<"patient ID of person rescreened and found -ve compared to "<<*iso_queue.begin()+q<<"ID number of person removed from Q"<<"\n";
//                                    }
//
//                                   iso_queue.erase(iso_queue.begin()+q); //dont know if this isactually erasing the right person!!!!!
//
//                                   if (DEFAULT_DEBUG==8){
//                                        cout<<"REMOVED FROM ISO_QUEUE >>>>>>>>>>>>>>>>>8"<<"\n";
//                                   }
//
//                                }
//                             }
//
//                //end of if they are screened correctly
//
//                } else{//not screened correctly
//
//                    if (DEFAULT_DEBUG==2){
//                        cout<<"believes C actually S"<<"  person="<<hos_pop[pick_patient].patientid<<"\n";
//                    }
//
//                    hos_pop[pick_patient].awareness_state = BELIEVED_C_ACTUALLY_S;
//
//                    //SCHEDULE THEIR RE-SCREEN
//                    future_intervention_events[time+rescreen_interval].push_back(newevent1);
//                    //ACT ON SCREEN RESULT
//                    future_intervention_events[time+TAT].push_back(newevent2);
//
//                            if (DEFAULT_DEBUG==2){
//                                cout<<"Act on screen scheduled for person="<<hos_pop[pick_patient].patientid<<"\n";
//                            }
//
//
//                            /**  (WRONGLY) CALLING INTERVENTION ROUTINES FROM SCREEN
//                            //PUT THEM IN ISOLATION
//                            isolation(pick_patient, rng);
//                            if (time>100){
//                            decolonization(rng, pick_patient);
//                            }
//                            */
//                        }//end of else - when they are not screened correctly
//
//            }
//
//
//            /**
//            else{
//                cout<<"disease state error, disease state is"<<int(hos_pop[pick_patient].disease_state)<<"\n";
//            }
//            */
//
//        //keeping tabs on number of screens taken




        ++no_ad_screens; //note that if more than one tye of screen is taken from the same patient they are counted sspearately here
        //cout<<no_ad_screens<<" number of screens count from admission screens"<<"\n";




        //  }//if screened
        } //end for(vector<screeningpolicycomponent>::...
}
    } //end of admission screening
//        else {//if not screened
//
//            if (DEFAULT_DEBUG==2){
//                cout<<"unknown"<<"  person="<<hos_pop[pick_patient].patientid<<"\n";
//            }
//
//            hos_pop[pick_patient].awareness_state = UNKNOWN;
//        }
//
//    }//end of screening
//

    void patients::scheduled_screening( all_patients* patient_ptr, gsl_rng *rng){
    //performs the scheduled screen for the specified patient (providing the individual has not yet been discharged form the ICU)
    //this schedules an ACT_ON_POS screen event if it is a positive screen and an ACT_ON_NEG screen event if it is a negative screen
        //cout<<"doing scheduled screening"<<"\n";
if (time>policy.timesteptoimplement){

        if (DEFAULT_DEBUG==4){
            cout <<"some scheduled screening is going on..........................."<<"\n";
        }


        struct events newevent1;
        struct events newevent2;
        struct events newevent3;

        newevent1.event_type=SCHEDULED_SCREEN;
        newevent2.event_type=ACT_ON_POS_SCREEN;
        newevent3.event_type=ACT_ON_NEG_SCREEN;
        newevent1.patient_ptr = patient_ptr;
        newevent2.patient_ptr = patient_ptr;
        newevent3.patient_ptr = patient_ptr;
        int rescreen_interval=7*timestepsperday; //specified in case  weekly rescreening is implemented

        int tat;           //turnaround time
        float sens, spec; //sensitivity and secificity

        float rannum, rannum6;
        bool screenthispatient; //patients are screened weekly with prob policy1.proportionscreenedweekly

        //set default values to those for conventional cultures
        //sens=conventional_culture.sensitivity;
        //spec=conventional_culture.specificity;
        //tat=conventional_culture.tat;

        if (DEFAULT_DEBUG==4){
            cout<<int (patient_ptr->disease_state)<<" actual disease state"<<"\n";
        }
        rannum = gsl_rng_uniform(rng);
        if(rannum<policy.proportionscreenedweekly) screenthispatient=true; else screenthispatient=false;//check to see if we select this patient
        //if not already discharged form ICU
        if (patient_ptr->ICU_no !=99 && screenthispatient){
            //schedule a screen for one weeks time

            future_intervention_events[time+rescreen_interval].push_back(newevent1);

            //then loop through components of the scheduled weekly screen
            for(vector<screeningpolicycomponent>::const_iterator iter=policy.weeklypostadmission.begin(); iter!=policy.weeklypostadmission.end();++iter){

                //keeping tabs on number of screens taken
                ++no_wkly_screens;


                //cout<<no_wkly_screens<<" number of screens count from weekly screens"<<"\n";
                if(iter->technology==CC) {sens=conventional_culture.sensitivity;spec=conventional_culture.specificity;tat=conventional_culture.tat;}
                if(iter->technology==CA) {sens=chromagar.sensitivity;spec=chromagar.specificity;tat=chromagar.tat;}
                if(iter->technology==CA_early) {sens=chromagar_early.sensitivity;spec=chromagar_early.specificity;tat=chromagar_early.tat;}
                if(iter->technology==PCR) {sens=pcr.sensitivity;spec=pcr.specificity;tat=pcr.tat;}
                if(iter->technology==I) {sens=ideal.sensitivity;spec=ideal.specificity;tat=ideal.tat;}

                if (DEFAULT_DEBUG==4){
                    cout<<patient_ptr->ICU_no <<"ICU number of re-screened patient"<<"\n";
                    cout<<"ID number of re-screened patient: " <<patient_ptr->patientid<<"\n";
                }




                if (patient_ptr->disease_state==COLONIZED||patient_ptr->disease_state==INFECTED){

                    rannum6 = gsl_rng_uniform(rng);


                    if (rannum6<sens){
                        if (DEFAULT_DEBUG==4){
                            cout<<"SCHEDULED SCREEN FOUND TRUE POSITIVe"<<"\n";
                        }

                        ++no_pos_screens;

                        if(iter->technology==CC){
                            ++no_pos_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_pos_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_pos_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }



                        future_intervention_events[time+tat].push_back(newevent2);//schedule acting on positive screen
                    } else {
                       if (DEFAULT_DEBUG==4){
                            cout<<"SCHEDULED SCREEN FOUND FALSE NEGATIVe"<<"\n";
                        }

                        ++no_neg_screens;

                        if(iter->technology==CC){
                            ++no_neg_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_neg_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_neg_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }



                        future_intervention_events[time+tat].push_back(newevent3);//schedule acting on negative screen
                    }
                } else if (patient_ptr->disease_state==SUSCEPTIBLE){
                    rannum6 = gsl_rng_uniform(rng);
                    if (rannum6<spec){
                        if (DEFAULT_DEBUG==4){
                            cout<<"SCHEDULED SCREEN FOUND TRUE NEGATIVE"<<"\n";
                        }

                        ++no_neg_screens;

                        if(iter->technology==CC){
                            ++no_neg_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_neg_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_neg_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }



                        future_intervention_events[time+tat].push_back(newevent3);//schedule acting on negative screen
                    } else {
                       if (DEFAULT_DEBUG==4){
                            cout<<"SCHEDULED SCREEN FOUND FALSE POSITIVe"<<"\n";
                        }

                        ++no_pos_screens;

                        if(iter->technology==CC){
                            ++no_pos_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_pos_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_pos_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }



                        future_intervention_events[time+tat].push_back(newevent2);//schedule acting on positive screen
                    }
                } else{  //if patient is not colonized or infected or susceptible  232 `
                    cout<<"disease_state not SUSCEPTIBLE, INFECTED or COLONIZED";
                    exit(0);
                }
            } //end for


        } //**MOBILE CODE** patient in ICU

}
    } //end patients::scheduled_screening

void patients::scheduled_dayofweek_screening( all_patients* patient_ptr, gsl_rng *rng){
    //performs the scheduled screen for the specified patient (providing the individual has not yet been discharged from the ICU)
    //this schedules an ACT_ON_POS screen event if it is a positive screen and an ACT_ON_NEG screen event if it is a negative screen
    //but unlike scheduled screening it does not schedule a new screen. Instead this routine implmenets screening policies that screen patients
    //on a particular day of the week. The procedure itself looks up to see if any regular screening is scheduled for the current
    //day of the week, and if implements it (if there is more than one type of screening for current day of week, all will
    //be implemented.
    //Note that the routine only implements the screening if the timestep is the first timestep of a new day. This is because
    // there may be more than one timestep per day, and under such circumstances we don't want to screen at every timestep
    //of the day

        //no newevent1 as that as to schedule a new screening event which we are not doing here

     if (time>policy.timesteptoimplement){


        struct events newevent2;
        struct events newevent3;

        newevent2.event_type=ACT_ON_POS_SCREEN;
        newevent3.event_type=ACT_ON_NEG_SCREEN;

        newevent2.patient_ptr = patient_ptr;
        newevent3.patient_ptr = patient_ptr;

        int currentday;

        int tat;           //turnaround time
        float sens, spec; //sensitivity and secificity

        float rannum, rannum6;
        bool screenthispatient; //patients are screened weekly with prob policy1.proportionscreenedweekly

        //set default values to those for conventional cultures
        //sens=conventional_culture.sensitivity;
        //spec=conventional_culture.specificity;
        //tat=conventional_culture.tat;

        currentday=get_weekday_number();
        rannum = gsl_rng_uniform(rng);
        if(rannum<policy.proportionscreenedweekly) screenthispatient=true; else screenthispatient=false;//check to see if we select this patient
        //if not already discharged from ICU and if it is the first time step of the day (since there could be >=1 time steps in one day)
        if (patient_ptr->ICU_no !=99 && is_it_first_time_step_of_day()==true && screenthispatient){

            //then loop through components of the scheduled weekday screens for the current day

            for(vector<screeningpolicycomponent>::const_iterator iter=policy.weekdayscreening[currentday].begin(); iter!=policy.weekdayscreening[currentday].end();++iter){

                //keeping tabs on number of screens taken
                ++no_wkly_screens;


                //cout<<no_wkly_screens<<" number of screens count from weekly screens"<<"\n";
                if(iter->technology==CC) {sens=conventional_culture.sensitivity;spec=conventional_culture.specificity;tat=conventional_culture.tat;}
                if(iter->technology==CA) {sens=chromagar.sensitivity;spec=chromagar.specificity;tat=chromagar.tat;}
                if(iter->technology==CA_early) {sens=chromagar_early.sensitivity;spec=chromagar_early.specificity;tat=chromagar_early.tat;}
                if(iter->technology==PCR) {sens=pcr.sensitivity;spec=pcr.specificity;tat=pcr.tat;}
                if(iter->technology==I) {sens=ideal.sensitivity;spec=ideal.specificity;tat=ideal.tat;}
                if (DEFAULT_DEBUG==4){
                    cout<<patient_ptr->ICU_no <<"ICU number of re-screened patient"<<"\n";
                    cout<<"ID number of re-screened patient: " <<patient_ptr->patientid<<"\n";
                }


                if (patient_ptr->disease_state==COLONIZED||patient_ptr->disease_state==INFECTED){

                    rannum6 = gsl_rng_uniform(rng);


                    if (rannum6<sens){
                        if (DEFAULT_DEBUG==4){
                            cout<<"SCHEDULED SCREEN FOUND TRUE POSITIVe"<<"\n";
                        }

                        ++no_pos_screens;

                        if(iter->technology==CC){
                            ++no_pos_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_pos_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_pos_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }


                        future_intervention_events[time+tat].push_back(newevent2);//schedule acting on positive screen
                    } else {

                       if (DEFAULT_DEBUG==4){
                            cout<<"SCHEDULED SCREEN FOUND FALSE NEGATIVe"<<"\n";
                        }

                        ++no_neg_screens;

                        if(iter->technology==CC){
                            ++no_neg_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_neg_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_neg_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }



                        future_intervention_events[time+tat].push_back(newevent3);//schedule acting on negative screen
                    }
                } else if (patient_ptr->disease_state==SUSCEPTIBLE){
                    rannum6 = gsl_rng_uniform(rng);
                    if (rannum6<spec){
                        if (DEFAULT_DEBUG==4){
                            cout<<"SCHEDULED SCREEN FOUND TRUE NEGATIVE"<<"\n";
                        }

                        ++no_neg_screens;

                        if(iter->technology==CC){
                            ++no_neg_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_neg_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_neg_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }



                        future_intervention_events[time+tat].push_back(newevent3);//schedule acting on negative screen
                    } else {
                       if (DEFAULT_DEBUG==4){
                            cout<<"SCHEDULED SCREEN FOUND FALSE POSITIVe"<<"\n";
                        }

                        ++no_pos_screens;

                        if(iter->technology==CC){
                            ++no_pos_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_pos_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_pos_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }



                        future_intervention_events[time+tat].push_back(newevent2);//schedule acting on positive screen
                    }
                } else{  //if patient is not colonized or infected or susceptible  232 `
                    cout<<"disease_state not SUSCEPTIBLE, INFECTED or COLONIZED";
                    exit(0);
                }
            } //end for


        } //**MOBILE CODE** patient in ICU and first time step of a new day

     }

    } //end patients::scheduled_dayofweek_screening




    void patients::clinical_screening( all_patients* patient_ptr, gsl_rng *rng){
    //performs the clinical for the specified infected patient (providing the individual has not yet been discharged form the ICU)
    //this schedules an ACT_ON_POS screen event if it is a positive screen and an ACT_ON_NEG screen event if it is a negative screen

if (time>policy.timesteptoimplement){

        if (DEFAULT_DEBUG==4){
            cout <<"some clinical screening is going on..........................."<<"\n";
        }


        struct events newevent1;
        struct events newevent2;

        newevent1.event_type=ACT_ON_POS_SCREEN;
        newevent2.event_type=ACT_ON_NEG_SCREEN;

        newevent1.patient_ptr = patient_ptr;
        newevent2.patient_ptr = patient_ptr;


        int tat;           //turnaround time
        float sens, spec; //sensitivity and secificity

        float rannum6;


        if (DEFAULT_DEBUG==4){
            cout<<int (patient_ptr->disease_state)<<" actual disease state"<<"\n";
        }

        //if not already discharged form ICU
        if (patient_ptr->ICU_no !=99){


            //then loop through components of the clinical screen
            for(vector<screeningpolicycomponent>::const_iterator iter=policy.clinicalscreening.begin(); iter!=policy.clinicalscreening.end();++iter){

                //keeping tabs on number of clinical screens taken
                ++no_clin_screens;


                //cout<<no_clin_screens<<" number of screens count from clinical screens"<<"\n";
                if(iter->technology==CC) {sens=conventional_culture.sensitivity;spec=conventional_culture.specificity;tat=conventional_culture.tat;}
                if(iter->technology==CA) {sens=chromagar.sensitivity;spec=chromagar.specificity;tat=chromagar.tat;}
                if(iter->technology==CA_early) {sens=chromagar_early.sensitivity;spec=chromagar_early.specificity;tat=chromagar_early.tat;}
                if(iter->technology==PCR) {sens=pcr.sensitivity;spec=pcr.specificity;tat=pcr.tat;}
                if(iter->technology==I) {sens=ideal.sensitivity;spec=ideal.specificity;tat=ideal.tat;}


                if (DEFAULT_DEBUG==4){
                    cout<<patient_ptr->ICU_no <<"ICU number of clinically screened patient"<<"\n";
                    cout<<"ID number of clinically screened patient: " <<patient_ptr->patientid<<"\n";
                }




                if (patient_ptr->disease_state==COLONIZED||patient_ptr->disease_state==INFECTED){

                    rannum6 = gsl_rng_uniform(rng);


                    if (rannum6<sens){
                        if (DEFAULT_DEBUG==4){
                            cout<<"SCHEDULED SCREEN FOUND TRUE POSITIVE"<<"\n";
                        }

                        ++no_pos_screens;

                        if(iter->technology==CC){
                            ++no_pos_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_pos_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_pos_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }

                        future_intervention_events[time+tat].push_back(newevent1);//schedule acting on positive screen
                    } else {
                       if (DEFAULT_DEBUG==4){
                            cout<<"SCHEDULED SCREEN FOUND FALSE NEGATIVE"<<"\n";
                        }

                        ++no_neg_screens;

                        if(iter->technology==CC){
                            ++no_neg_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_neg_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_neg_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }



                        future_intervention_events[time+tat].push_back(newevent2);//schedule acting on negative screen
                    }
                } else if (patient_ptr->disease_state==SUSCEPTIBLE){
                    rannum6 = gsl_rng_uniform(rng);
                    if (rannum6<spec){
                        if (DEFAULT_DEBUG==4){
                            cout<<"SCHEDULED SCREEN FOUND TRUE NEGATIVE"<<"\n";
                        }

                        ++no_neg_screens;

                        if(iter->technology==CC){
                            ++no_neg_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_neg_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_neg_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }



                        future_intervention_events[time+tat].push_back(newevent2);//schedule acting on negative screen
                    } else {
                       if (DEFAULT_DEBUG==4){
                            cout<<"SCHEDULED SCREEN FOUND FALSE POSITIVE"<<"\n";
                        }

                        ++no_pos_screens;

                        if(iter->technology==CC){
                            ++no_pos_CC_screens;
                        }
                        if(iter->technology==CA){
                            ++no_pos_CA_screens;
                        }
                        if(iter->technology==CA_early){
                            ++no_pos_CA_early_screens;
                        }
                        if(iter->technology==PCR){
                            ++no_PCR_screens;
                        }
                        if(iter->technology==I){
                            ++no_ideal_screens;
                        }



                        future_intervention_events[time+tat].push_back(newevent1);//schedule acting on positive screen
                    }
                } else{  //if patient is not colonized or infected or susceptible  232 `
                    cout<<"disease_state not SUSCEPTIBLE, INFECTED or COLONIZED";
                    exit(0);
                }
            } //end for


        } //**MOBILE CODE** patient in ICU

}

    } //end patients::clinical_screening



    void patients::remove_control_measures_for_those_believed_negative(all_patients* patient_ptr, gsl_rng *rng){

// cout<<"removed control measures for person who is " <<int(patient_ptr->disease_state)<<"on day "<< time<<"\n";


//remove isolation
        if (patient_ptr->isolation==true){ //i.e. only do this if currently isolated
            patient_ptr->isolation=false;
            --ISO;

//cout<<"patient "<<patient_ptr->patientid << " is being taken out of isolation in light of new results on day " <<time <<"\n";

            patient_ptr->IQ= IQ ;//baseline IQ
        }


//remove secondary isolation
        if (patient_ptr->secndisolation==true){ //i.e. only do this if currently isolated
            patient_ptr->secndisolation=false;
            --SECISO;
            patient_ptr->IQ= IQ ;//baseline IQ
        }


        //if a person in the iso_queue is believed negative take them off queue and stop secondary isolation for them
        for (unsigned long int q=0; q<iso_queue.size(); q++){
            if (iso_queue[q]==patient_ptr->patientid){

                iso_queue.erase(iso_queue.begin()+q);
                patient_ptr->secndisolation=false;
                patient_ptr->IQ=IQ;
                --SECISO;

                if (DEFAULT_DEBUG==8){
                    cout<<*iso_queue.begin()+q<<"person taken off iso Q because they were belived -ve"<<"\n";
                }

            }
         }


            //if removals have created space in the IW
            //and if there are people waiting (i.e. in iso_queue)
            //move them in (and stop their secondary isolation)

        if (ISO<policy.primary_isolation_capacity && iso_queue.size()>0){


            if (DEFAULT_DEBUG==8){
               cout<<*iso_queue.begin()<<"person being moved from Q into ISO"<<"\n";
            }

            hos_pop[*iso_queue.begin()].isolation=true;
            ISO++;
            patient_ptr->secndisolation=false;
            --SECISO;

            hos_pop[*iso_queue.begin()].IQ = IQ * (1-policy.primary_isolation_effectiveness );


            iso_queue.erase(iso_queue.begin());

            if (DEFAULT_DEBUG==8){
                cout<<"REMOVED FROM ISO_QUEUE >>>>>>>>>>>>>>>>>2"<<"\n";
            }
        }


//stop decolonization
        if (patient_ptr->decolonization==true){ //i.e. only do this if the patient is in the middle of undergoing treatment when they are found to be negative
			end_decolonization(patient_ptr, rng);
        }







    } //end of remove control measures for those believed to be negative

    void patients::implement_control_measures_for_positives(all_patients* patient_ptr, gsl_rng *rng){



    //implements whatever control measures are in place for the selected patient under the assumption that they
     //are believed to be MRSA+
    //interventions can result in changes in the infectivity (IQ) of those patients either known or believed to be +ve through screening

    //currently this just does patient isolation...need to add decolonization

        if(policy.isolateifpositive==true){ //isolate of policy is to isolate positives and if there is capacity

            //if there is space
            if (ISO<policy.primary_isolation_capacity){
                            if (patient_ptr->isolation==false){ //i.e. only do this if not currently isolated

                                //cout<<"patient "<<patient_ptr->patientid << " is being isolated when they werent before on day " <<time <<"\n";

                                patient_ptr->isolation=true;
                                ISO++;
                                patient_ptr->IQ= IQ * (1- policy.primary_isolation_effectiveness);


                                if (DEFAULT_DEBUG==8){
                                    cout<<patient_ptr->patientid<< "  ID no of isolated patient"<<"\n";
                                }
                            }
            }  else {        //if there is no space, put them on the queue and put them under secondary isolation

                int alreadyinisoqueuecount = 0;

                for (unsigned long int q=0; q<iso_queue.size(); q++){
                    if (patient_ptr->patientid==iso_queue[q]){
                        alreadyinisoqueuecount++;
                    }
                }

                if (alreadyinisoqueuecount>0){
                        if (DEFAULT_DEBUG==8){
                            cout<<"already in isoqueue"<<"\n";
                        }
                } else{ //not already in the queue

                    if (DEFAULT_DEBUG==8){
                        cout<<iso_queue.size() <<"size of iso_queue before push_back"<<"\n";
                    }

                    iso_queue.push_back(patient_ptr->patientid);

                    patient_ptr->secndisolation=true;
                    ++SECISO;
                    patient_ptr->IQ=IQ* (1-policy.secondary_isolation_effectiveness);




                        if (DEFAULT_DEBUG==8){
                            cout<<patient_ptr->patientid<< "  ID no of patient who was to be isolated but instead put on queue"<<"\n";
                            cout<<iso_queue.size() <<"size of iso_queue after push_back"<<"\n";
                            for (unsigned long int i = 0; i<iso_queue.size(); ++i){
                                cout<<iso_queue[i] <<"ID number of patient in isoqueue"<<"\n";
                            }
                        }
                } //**MOBILE CODE** (alreadyinisoqueue)
            }//**MOBILE CODE** (ISO<Isolation Capacity)
        }//**MOBILE CODE** (policy.isolateifpositive)


        if(policy.decolifpositive==true){
			if(patient_ptr->decolonization==false){
                decol(patient_ptr, rng);
                ++decolcount;
			}

        }


    }//end implement_control_measures_for_positives




    void patients::implement_control_measures_for_all(gsl_rng *rng, int pick_patient){

    //implements whatever control measures are in place for admitted patients

  if (time>policy.timesteptoimplement){

        if(policy.preemptiveisolationforall==true){


           //if there is space
            if (ISO<policy.primary_isolation_capacity){
                            if (hos_pop[pick_patient].isolation==false){ //i.e. only do this if not currently isolated
                                hos_pop[pick_patient].isolation=true;
                                ISO++;

                                //cout<<"patient "<<hos_pop[pick_patient].patientid << " is being isolated (due to control measures for all) when they werent before on day " <<time <<"\n";

                                hos_pop[pick_patient].IQ= IQ * (1-policy.primary_isolation_effectiveness);

                                if (DEFAULT_DEBUG==8){
                                    cout<<hos_pop[pick_patient].patientid<< "  ID no of isolated patient"<<"\n";
                                }


                            }
            }  else {        //if there is no space, put them on the queue and pu them under secondary isolation

                int alreadyinisoqueuecount = 0;

                for (unsigned long int q=0; q<iso_queue.size(); q++){
                    if (hos_pop[pick_patient].patientid==iso_queue[q]){
                        alreadyinisoqueuecount++;
                    }
                }

                if (alreadyinisoqueuecount>0){
                        if (DEFAULT_DEBUG==8){
                            cout<<"already in isoqueue"<<"\n";
                        }
                } else{ //not already in the queue

                    if (DEFAULT_DEBUG==8){
                        cout<<iso_queue.size() <<"size of iso_queue before push_back"<<"\n";
                    }

                    iso_queue.push_back(hos_pop[pick_patient].patientid);

                    hos_pop[pick_patient].secndisolation=true;
                    ++SECISO;
                    hos_pop[pick_patient].IQ=IQ* (1-policy.secondary_isolation_effectiveness );


                        if (DEFAULT_DEBUG==8){
                            cout<<hos_pop[pick_patient].patientid<< "  ID no of patient who was to be isolated but instead put on queue"<<"\n";
                            cout<<iso_queue.size() <<"size of iso_queue after push_back"<<"\n";
                            for (unsigned long int i = 0; i<iso_queue.size(); ++i){
                                cout<<iso_queue[i] <<"ID number of patient in isoqueue"<<"\n";
                            }
                        }
                } //**MOBILE CODE** (alreadyinisoqueue)
            }//**MOBILE CODE** (ISO<Isolation Capacity)


        }//end     if(policy.preemptiveisolationforall){


        if(policy.blanketdecolforall==true){
			if (hos_pop[pick_patient].decolonization==false){
                decol(rng, pick_patient);
                ++decolcount;
			}
        }
        }
    }//end of implement control measures for all




    void patients::implement_control_measures_for_high_risk(gsl_rng *rng, int pick_patient){

    //implements whatever control measures are in place for high risk admitted patients



      if (time>policy.timesteptoimplement){
        if(hos_pop[pick_patient].risk_group==HIGH){






            if(policy.preemptiveisolationforhighrisk==true){


           //if there is space
            if (ISO<policy.primary_isolation_capacity){
                            if (hos_pop[pick_patient].isolation==false){ //i.e. only do this if not currently isolated
                                hos_pop[pick_patient].isolation=true;
                                ISO++;
                                hos_pop[pick_patient].IQ= 1 * (1-policy.primary_isolation_effectiveness);


                                if (DEFAULT_DEBUG==8){
                                    cout<<hos_pop[pick_patient].patientid<< "  ID no of isolated patient"<<"\n";
                                }
                            }
            }  else {        //if there is no space, put them on the queue and in mean time put them under secondary isolation

                int alreadyinisoqueuecount = 0;

                for (unsigned long int q=0; q<iso_queue.size(); q++){
                    if (hos_pop[pick_patient].patientid==iso_queue[q]){
                        alreadyinisoqueuecount++;
                    }
                }

                if (alreadyinisoqueuecount>0){
                        if (DEFAULT_DEBUG==8){
                            cout<<"already in isoqueue"<<"\n";
                        }
                } else{ //not already in the queue

                    if (DEFAULT_DEBUG==8){
                        cout<<iso_queue.size() <<"size of iso_queue before push_back"<<"\n";
                    }

                    iso_queue.push_back(hos_pop[pick_patient].patientid);


                    hos_pop[pick_patient].secndisolation=true;
                    ++SECISO;
                    hos_pop[pick_patient].IQ=IQ* (1-policy.secondary_isolation_effectiveness);


                        if (DEFAULT_DEBUG==8){
                            cout<<hos_pop[pick_patient].patientid<< "  ID no of patient who was to be isolated but instead put on queue"<<"\n";
                            cout<<iso_queue.size() <<"size of iso_queue after push_back"<<"\n";
                            for (unsigned long int i = 0; i<iso_queue.size(); ++i){
                                cout<<iso_queue[i] <<"ID number of patient in isoqueue"<<"\n";
                            }
                        }
                } //**MOBILE CODE** (alreadyinisoqueue)
            }//**MOBILE CODE** (ISO<Isolation Capacity)


        }//end     if(policy.preemptiveisolationforall){



        if(policy.decolforhighrisk==true){
			if (hos_pop[pick_patient].decolonization==false){ //need to check here a we only increment decolcount if not currently an therapy
                decol(rng, pick_patient);
                ++decolcount;
			}
        }

        }


        }
    }//end of implement control measures for high risk





    void patients::decol(all_patients* patient_ptr, gsl_rng *rng){


//called from screening routine and therefore using patient pointer

     // to implement need to set decolonization_mup or decolonization_mup_chx
     // and need to schedule end of treatment.
     // Note that duration of treatment & effect of treatment on susceptibility and
     // transmissibility should not be specified in this routine, but in the same place
     // as all the other parameters - and value of duration of treatment and impact on
     // transmission parameters should be customizable at run-time (instead of only at comile time as currently

        if (DEFAULT_DEBUG==64){
            cout<<"decolonization routine called"<<"\n";
        }


        struct events newevent;
        newevent.event_type=FINISH_TREATMENT;
        newevent.patient_ptr = patient_ptr;

		if (patient_ptr->decolonization==false){
         if(patient_ptr->disease_state==SUSCEPTIBLE){
             patient_ptr->Pdc=patient_ptr->Pdc*(1-policy.bodywash_pdc_effectiveness);
             patient_ptr->Pdi=patient_ptr->Pdi*(1-policy.bodywash_pdi_effectiveness);
             patient_ptr->IQ=patient_ptr->IQ*(1-policy.bodywash_IQ_effectiveness);
             patient_ptr->ProgProb=patient_ptr->ProgProb*(1-policy.bodywash_progprob_effectiveness);
         }



         if(patient_ptr->disease_state==COLONIZED){
			 patient_ptr->Pdc=patient_ptr->Pdc*(1-policy.bodywash_pdc_effectiveness);
             patient_ptr->IQ=patient_ptr->IQ*(1-policy.bodywash_IQ_effectiveness);
             patient_ptr->ProgProb=patient_ptr->ProgProb*(1-policy.bodywash_progprob_effectiveness);
             patient_ptr->Pdi=patient_ptr->Pdi*(1-policy.bodywash_pdi_effectiveness);

         }
		 patient_ptr->decolonization=true;

		}
        future_intervention_events[time+length_of_treatment].push_back(newevent);//schedule end of treatment



    }



     void patients::process_readmission(gsl_rng * rng){


//        unsigned long int pick_patient;
        float rannum;
        rannum = gsl_rng_uniform(rng);//uniform random number: includes 0 but excludes 1



        /* int ICUsize;

         for (int m=0; m<DEFAULT_ICU_SIZE; m++){


                    if (ICU_patients[m].ICU_patient_ptr!=NULL){

                        ++ICUsize;

                    }
         }


*/
        // if (ICUsize<DEFAULT_ICU_SIZE){


         if (DEFAULT_DEBUG==72){
                cout<<*admis_queue.begin()<<" Person being moved from admission Q into ICU"<<"\n";
                cout<<admis_queue.size()<<" Size of admission Q"<<"\n";
         }



                    for (int m=0; m<DEFAULT_ICU_SIZE; m++){


                    if (ICU_patients[m].ICU_patient_ptr==NULL){           ///what about if there is more than- 1 empty bed?

                            if (DEFAULT_DEBUG==72){
                                //cout<<ICU_patients[m].bed_no<<"bed_no_readmission";
                            }

                            ICU_patients[m].ICU_patient_ptr=&hos_pop[*admis_queue.begin()];

                            ICU_patients[m].ICU_patient_ptr->ICU_no=m;

                            ICU_patients[m].ICU_patient_ptr->specialty=ICU;

                            ICU_patients[m].ICU_patient_ptr->risk_group=HIGH;

                            ICU_patients[m].ICU_patient_ptr->no_days_in_hos=0;

                            if (ICU_patients[m].ICU_patient_ptr->disease_state==COLONIZED)

                            {

                                --S;
                                ++C;

                                }

                              if (ICU_patients[m].ICU_patient_ptr->disease_state==INFECTED)

                            {

                                --S;
                                ++I;

                                }


                                   if (DEFAULT_DEBUG==72){
                                cout<<ICU_patients[m].ICU_patient_ptr->ICU_no<<"  ICU_no after re-admission"<<"\n";
                                cout<<int(ICU_patients[m].ICU_patient_ptr->disease_state)<<"  disease state of re-admittee"<<"\n";
                                cout<<ICU_patients[m].ICU_patient_ptr->no_hos<<"  number of hospitalisations (re-admission)"<<"\n";
                            }


                    }//if bed is empty

                            }//for whole ICU


            admis_queue.erase (admis_queue.begin());//if someone is admitted from the admission queue remove them from the queue




            //}


            /* else {        //if there is no space, put them on the queue

            int alreadyinadmisqueuecount = 0;

              for (unsigned long int q=0; q<admis_queue.size(); q++){
                    if (patient_ptr->patientid==admis_queue[q]){
                       alreadyinadmisqueuecount++;
                    }
                }



                if (alreadyinadmisqueuecount>0){
                        if (DEFAULT_DEBUG==8){
                            cout<<"already in isoqueue"<<"\n";
                        }
                } else{ //not already in the queue

                    if (DEFAULT_DEBUG==8){
                        cout<<admis_queue.size() <<"size of iso_queue before push_back"<<"\n";
                    }

                    admis_queue.push_back(patient_ptr->patientid);
                }

            }*/

         }








    void patients::end_decolonization( all_patients* patient_ptr, gsl_rng *rng){


        //end decolonization therapy, reset patient's IQ etc, and decolonize patient with some probability
		if(patient_ptr->decolonization==true){
			patient_ptr->IQ=	patient_ptr->baselineIQ;
			patient_ptr->Pdc=patient_ptr->baselinePdc;
			patient_ptr->Pdi=patient_ptr->baselinePdi;
			patient_ptr->ProgProb=patient_ptr->baselineProgProb;


//		 patient_ptr->IQ=IQ; //to get baseline IQ from dist: get_IQ(rng);
//         patient_ptr->Pdc=patient_ptr->Pdc/(1-policy.bodywash_pdc_effectiveness);
//         patient_ptr->Pdi=patient_ptr->Pdi/(1-policy.bodywash_pdi_effectiveness);
//         patient_ptr->ProgProb=patient_ptr->ProgProb/(1-policy.bodywash_progprob_effectiveness);


        //assume at the end of the treatment period a certain proportion of individuals have successful treatment i.e. become S
		 patient_ptr->decolonization=false;

        float random;
        random = gsl_rng_uniform(rng);

                if (random<PropTreatmentSuccessful){
                    if (patient_ptr->disease_state==COLONIZED){
                        patient_ptr->disease_state = SUSCEPTIBLE;
                        --C;
                        ++S;
                    }
                }

	   }
    }// end of end decolonization






    void patients::decol(gsl_rng *rng, int pick_patient){


 //called from admission routine therefore using pick patient

     // to implement need to set decolonization_mup or decolonization_mup_chx
     // and need to schedule end of treatment.
     // Note that duration of treatment & effect of treatment on susceptibility and
     // transmissibility should not be specified in this routine, but in the same place
     // as all the other parameters - and value of duration of treatment and impact on
     // transmission parameters should be customizable at run-time

        if (DEFAULT_DEBUG==64){
            cout<<"decolonization routine called"<<"\n";
        }

        struct events newevent;
        newevent.event_type=FINISH_TREATMENT;
        newevent.patient_ptr = &hos_pop[pick_patient];
        //int length_of_treatment=5;//treatment lasts 5 days then stopped
        //float DecolEffectonPdc = 0.72;// Barry's guess
        //float DecolEffectonIQ =0.72;//from Barry's guess
        //float DecolEffectonProgProb =0.5;//made up

		if (hos_pop[pick_patient].decolonization==false){

         if(hos_pop[pick_patient].disease_state==SUSCEPTIBLE){
             hos_pop[pick_patient].IQ=hos_pop[pick_patient].IQ*(1-policy.bodywash_IQ_effectiveness);
             hos_pop[pick_patient].ProgProb=hos_pop[pick_patient].ProgProb*(1-policy.bodywash_progprob_effectiveness);//*0.5);//made up
             hos_pop[pick_patient].Pdc=hos_pop[pick_patient].Pdc*(1-policy.bodywash_pdc_effectiveness);
             hos_pop[pick_patient].Pdi=hos_pop[pick_patient].Pdi*(1-policy.bodywash_pdi_effectiveness);
         }

         if(hos_pop[pick_patient].disease_state==COLONIZED){
             hos_pop[pick_patient].IQ=hos_pop[pick_patient].IQ*(1-policy.bodywash_IQ_effectiveness);
             hos_pop[pick_patient].ProgProb=hos_pop[pick_patient].ProgProb*(1-policy.bodywash_progprob_effectiveness);//*0.5);//made up
             hos_pop[pick_patient].Pdi=hos_pop[pick_patient].Pdi*(1-policy.bodywash_pdi_effectiveness);
             hos_pop[pick_patient].Pdc=hos_pop[pick_patient].Pdc*(1-policy.bodywash_pdc_effectiveness);

         }
			hos_pop[pick_patient].decolonization=true;
         future_intervention_events[time+length_of_treatment].push_back(newevent);//schedule end of treatment

		}

    }


void patients::print_events(){

 char homeduration_file[]="homeduration.txt";
 std::ofstream durStream(homeduration_file);


for(std::map<int, unsigned long int>::const_iterator iter=duration_list.begin(); iter!=duration_list.end();++iter){

durStream<<(*iter).first<<" \t"<<(*iter).second<<" \n";

        }


  cout<<"\n Scheduled movement events : \n";

  for( vector<events>::const_iterator iter=future_movement_events[time].begin(); iter!=future_movement_events[time].end(); ++iter){

    cout<<int(iter->event_type)<<" = iter-> event type "<<"\n";
   }

  cout<<"\n Scheduled infection events : \n";

  for(vector<events>::const_iterator iter=future_infection_events[time].begin(); iter!=future_infection_events[time].end(); ++iter){
    cout<<int(iter->event_type)<<" = iter-> event type"  <<"\n";
  }

  cout<<"\n Scheduled intervention events : \n";
   for( vector<events>::const_iterator iter=future_intervention_events[time].begin(); iter!=future_intervention_events[time].end(); ++iter){
  cout<<int(iter->event_type)<<" = iter-> event type "<<"\n";

   }
}


void patients::perform_events(gsl_rng *rng){

        // performs all scheduled events and schedules new events

        //OUTPUTS DAILY COUNTS: no of S->C , S->I , C->I events
        //                      no of discharges, admissions and deaths

        //cout<<"perform events"<<"\n";

        int StoC=0;
        int StoI=0;
        int CtoI=0;
        int Conadmission=0;

        int discharges=0;
        int admissions=0;
        int readmissions=0;
        int deaths = 0;
//        int readmission = 0;


        struct events newevent;
        struct events newevent1;
      //  struct events newevent2;

        T_AwarenessState CorIstate1, CorIstate2;
        T_AwarenessState old_awareness_state;
    //****************************MOVEMENT EVENTS***********************************

          for( vector<events>::const_iterator iter=future_movement_events[time].begin(); iter!=future_movement_events[time].end(); ++iter){

                if (DEFAULT_DEBUG==72){
                    cout<<int(iter->event_type)<<" = iter-> event type "<<int(time)<<" time "<<"\n";
                }

                switch(int(iter->event_type)){

                        case ICU_DIS_AND_ADMIT:

                            if  (iter->patient_ptr->ICU_no!=99){

                                //cout<<"discharge patient ID"<<int(iter->patient_ptr->patientid)<<"\n";

                                if (DEFAULT_DEBUG==72){
                                    cout<<iter->patient_ptr<<"   address of patient to be discharged"<<"\n";
                                }

                                //cout<< "DOING A DISCHARGE to patient "<<iter->patient_ptr->ICU_no<<"\n";

                                cumulative_good_bd = cumulative_good_bd + (iter->patient_ptr->no_days_in_hos);


                                //cout<<iter->patient_ptr->no_days_in_hos<< " "<<"\n";

                                process_ICU_discharges(iter->patient_ptr, rng);//call process_ICU_discharges //added this rng because calling effect_of_ICU (from process_ICU discharges) required a rnadomnumber, so had to add it to function requirements of process_ICU_discharges
                                ++discharges;
                                ++cumulative_discharges;

                                 if (iter->patient_ptr->decolonization==true){ //i.e. only do this if the patient is in the middle of undergoing treatment when they die
	        //cout<<"xxxx ending decol for a discharge";
			end_decolonization(iter->patient_ptr, rng);
                                 }



                                if  (admis_queue.size()>0){

            if (DEFAULT_DEBUG==72){
                                    cout<<"   patient re-admitted"<<"\n";
            }



           process_readmission(rng);
           ++readmissions;
           ++cumulative_readmissions;

       if (hos_pop[*admis_queue.begin()].disease_state==COLONIZED){

               ++Conadmission;
               ++cumulative_Conadmission;

           }

            }

                                else{

                                     if (DEFAULT_DEBUG==72){
                                    cout<<"   patient admitted"<<"\n";
                                }



                               process_ICU_admissions(rng);
                               ++admissions;
                               ++cumulative_admissions;

                                }




                            }


                        break;


                        case DEATH:

                            if  (iter->patient_ptr->ICU_no!=99){

                                if (DEFAULT_DEBUG==1){
                                    cout<<iter->patient_ptr<<"   address of patient that dies"<<"\n";
                                }

                    //cout<< "DOING A DEATH to patient "<<iter->patient_ptr->ICU_no<<"\n";

                                process_ICU_deaths(iter->patient_ptr, rng);
                                ++deaths;
                                ++cumulative_deaths;

                                 if (int(iter->event_type) == READMISSION){


                               process_readmission( rng);
                               ++readmissions;
                               ++cumulative_readmissions;

                                }

                                else
                                {


                               process_ICU_admissions(rng);
                               ++admissions;
                               ++cumulative_admissions;

                                }

                            }


                             break;

                             if  (iter->patient_ptr->ICU_no=99){
                             case READMISSION:


                            admis_queue.push_back(iter->patient_ptr->patientid);

                      if (DEFAULT_DEBUG==72){

                    cout<<iter->patient_ptr<<"   address of patient to be re-admitted"<<"\n";

                    cout<<admis_queue.size() <<"size of admission queue after push_back"<<"\n";

                      }

                             }
                       break;


                              /*if (iter->patient_ptr->specialty = COM)
                               {

                                //cout<<"discharge patient ID"<<int(iter->patient_ptr->patientid)<<"\n";



                                //cout<< "DOING A DISCHARGE to patient "<<iter->patient_ptr->ICU_no<<"\n";

                              //  cumulative_good_bd = cumulative_good_bd + (iter->patient_ptr->no_days_in_hos);


                                //cout<<iter->patient_ptr->no_days_in_hos<< " "<<"\n";


                                //{
                                // if (DEFAULT_DEBUG==1){
                                  //  cout<<iter->patient_ptr<<"   address of patient that is readmitted"<<"\n";
                             //  }



                    //iter->patient_ptr->waitingforadmission=true;

                       if (DEFAULT_DEBUG==72){
                                    cout<<iter->patient_ptr<<"   address of patient to be re-admitted"<<"\n";
                                }


                             //    process_readmission(iter->patient_ptr, rng);
                             //    ++readmission;
                             //    ++cumulative_readmissions;

                            }*/


                }//end of switch

            }//end of for movements events

          //}//for movement events
    //****************************INFECTION EVENTS**********************************



        for(vector<events>::const_iterator iter=future_infection_events[time].begin(); iter!=future_infection_events[time].end(); ++iter){

            switch(int(iter->event_type)){

            case COLONIZE:
            //code to colonize (if susceptible) and schedule recovery



                if  (iter->patient_ptr->ICU_no!=99){ //i.e. if patient is in the ICU

                    if(iter->patient_ptr->disease_state==SUSCEPTIBLE){
                        iter->patient_ptr->disease_state=COLONIZED;
                        --S;
                        ++C;

                        ++StoC;
                        ++cumulative_StoC;

                        //cout<<iter->patient_ptr->patientid<< "id of patient becoming col **********" <<"\n";

                        iter->patient_ptr->timecolonized=time;

                        newevent.event_type=RECOVERFROMCOL;//scheduling
                        newevent.patient_ptr=iter->patient_ptr;

                        if (DEFAULT_DEBUG==32){
                            cout<<newevent.patient_ptr<<" address of person who has been added  onto list"<<"\n";
                        }

//@@@ col_duration paramets must also be sPecified centrally, and can calculate actual value using a function call.
// leaving it to the function to draw a random number from a distribution if required.

                        int col_duration ;//10;//4; //needs to be drawn from dist
                        col_duration =get_col_duration(rng);

                        future_infection_events[time+col_duration].push_back(newevent);

                        if (DEFAULT_DEBUG==32){
                         cout<< "colonization occurred to patient   "<<iter->patient_ptr->patientid<<"\n";
                        }

                    }//end of if susceptible
                }//end of if they are in the ICU


            break;


            case INFECT:
            //code to infect susceptible inidividuals and schedule their recovery



            if  (iter->patient_ptr->ICU_no!=99){



                if(iter->patient_ptr->disease_state==SUSCEPTIBLE){
                    iter->patient_ptr->disease_state=INFECTED;
                    iter->patient_ptr->timeinfected=time;
                    --S;
                    ++I;

                    ++StoI; ++cumulative_StoI;

                    //cout<<iter->patient_ptr->patientid<< "id of patient becoming inf *********" <<"\n";

                    iter->patient_ptr->timeinfected=time;

                    //queue up a recovery
                    newevent.event_type=RECOVERFROMINF;//scheduling
                    newevent.patient_ptr=iter->patient_ptr;

                    int inf_duration; //to be drawn from a distribution.
                    inf_duration =get_inf_duration(rng);
                    future_infection_events[time+inf_duration].push_back(newevent);



                    // queue up the clinical swab
                    newevent1.event_type=CLINICAL_SCREEN;
                    newevent1.patient_ptr=iter->patient_ptr;
                    future_intervention_events[time+policy.delaybeforeclinicalswab].push_back(newevent1);



                    if (DEFAULT_DEBUG==32){
                        cout<<"infection occured to patient   "<<iter->patient_ptr->patientid<< "\n";
                    }
                }//end of if they are sus



                if(iter->patient_ptr->disease_state==COLONIZED){
                    iter->patient_ptr->disease_state=INFECTED;
                    iter->patient_ptr->timeinfected=time;
                    --C;
                    ++I;

                    ++CtoI;
                    ++cumulative_CtoI;

                    iter->patient_ptr->timeinfected=time;


                    newevent.event_type=RECOVERFROMINF;//scheduling
                    newevent.patient_ptr=iter->patient_ptr;

                    int inf_duration; // to be drawn from a distribution.
                    inf_duration =get_inf_duration(rng);
                    future_infection_events[time+inf_duration].push_back(newevent);



                    // queue up the clinical swab
                    newevent1.event_type=CLINICAL_SCREEN;
                    newevent1.patient_ptr=iter->patient_ptr;
                    future_intervention_events[time+policy.delaybeforeclinicalswab].push_back(newevent1);



                    if (DEFAULT_DEBUG==32){
                        cout<<"transmission to an already colonized individual occurred to patient   "<<iter->patient_ptr->patientid<<"\n";
                    }
                }



            }//end of if they are in the ICU



            break;

            case RECOVERFROMCOL:
            //recovery of a colonized individual

            if  (iter->patient_ptr->ICU_no!=99){ //i.e. if patient is in the ICU

                if(iter->patient_ptr->disease_state==COLONIZED){
                    iter->patient_ptr->disease_state=SUSCEPTIBLE;
                    --C;
                    ++S;

                    iter->patient_ptr->timerecovered=time;
                    //iter->patient_ptr->extra_stay_to_be_added=0;

                    if (DEFAULT_DEBUG==32){
                      cout<<"recovery occured from col patient    "<<iter->patient_ptr->patientid<<"\n";
                    }

                }

                //if the colonized person had become inf in the meantime then recovery from col wont occur
                else

                    if (DEFAULT_DEBUG==32){
                      cout<< "patient sceduled to recover from colonized state is no longer in colonized state " <<"\n";
                    }
                }

            break;


            case RECOVERFROMINF:
                //recovery of an infected individual
                if  (iter->patient_ptr->ICU_no!=99){

                    if(iter->patient_ptr->disease_state==INFECTED){

                        iter->patient_ptr->disease_state=SUSCEPTIBLE;
                        --I;
                        ++S;

                        iter->patient_ptr->timerecovered=time;
                    //iter->patient_ptr->extra_stay_to_be_added=0;

                        if (DEFAULT_DEBUG==32){
                            cout<<"recovery occured from inf to patient  "<<iter->patient_ptr->patientid<<"\n";
                        }
                    }
                }

            break;


            default:
    //
                if (DEFAULT_DEBUG==32){
                     cout<<"Inappropriate event type found by perform_events. Person "<<iter->patient_ptr->patientid<<" & event type "<<int(iter->event_type)<<"\n";
                }


            break;

            }//switch

        }//for (infection events)


      future_infection_events.erase (time);//can delete events scheduled for current time as they have all now been dealt with




    //****************************INTERVENTION EVENTS**********************************

        for( vector<events>::const_iterator iter=future_intervention_events[time].begin(); iter!=future_intervention_events[time].end(); ++iter){


            switch(int(iter->event_type)){

            case ACT_ON_POS_SCREEN:
                //cout<<"\n *** reached case ACT_ON_POS_SCREEN ****** \n";
                //this should reset the awareness state and implement any control measures which are in place
                if  (iter->patient_ptr->ICU_no!=99 ){//check that they are in the ICU before doing anything to them - dont want to isolate someone not actually in the ward
                    //first change awareness state
                    if(iter->patient_ptr->disease_state==COLONIZED){
                      iter->patient_ptr->awareness_state=KNOWN_C;
                    } else if(iter->patient_ptr->disease_state==INFECTED){
                      iter->patient_ptr->awareness_state=KNOWN_I;
                    } else if(iter->patient_ptr->disease_state==SUSCEPTIBLE){
                        iter->patient_ptr->awareness_state=BELIEVED_C_ACTUALLY_S; //since either a false positve or subsequently cleared
                    } else {
                        cout<<"Undefined disease_state used";
                        exit(0);
                    }

                  //subsequent action depends on control policy in place
                        implement_control_measures_for_positives(iter->patient_ptr, rng); //deals with all control measurs for those believed mrsa+
                    // isolation(iter->patient_ptr, rng);
                    //decolonization(iter->patient_ptr, rng);

                }
                iter->patient_ptr->previous_screen_result=POSITIVE;
                iter->patient_ptr->num_consec_neg_screens_following_a_pos=0;
                iter->patient_ptr->everpositiveswabthisadmission=true;
                //hos_pop[i].num_consec_neg_screens_following_a_pos=0;
    //
    //    if (iter->patient_ptr->ICU_no==99){
    //            cout<<"TRYING TO ACT ON SCREEN OF SOMEONE WHO HAS ALREADY BEEN DISCHARGED"<<"\n";

    //            }



            break;

            case ACT_ON_NEG_SCREEN:

              old_awareness_state=iter->patient_ptr->awareness_state;
                //this might reset the awareness state (depending on how many negatives are required before confirmed neg
                // and control measures which are in place may change as a result of the changed awareness state
                //basic ideas is that if it is the first swab, or if no positive swabs have been received
                // during current stay believed status is negative. If there have been prior positive swabs during current stay
                // believed status is only negative if there have been a certain number of consecutive negative swabs
                // (this required num

                if(iter->patient_ptr->everpositiveswabthisadmission){//i.e. if they have at least one prior +ve culture this admission
                    ++iter->patient_ptr->num_consec_neg_screens_following_a_pos; //number of consecutive negative screens so far.
                }
                if  (iter->patient_ptr->ICU_no!=99 ){//check that they are in the ICU before doing anything to them - dont want to isolate someone not actually in the ward
                    //first change awareness state
                    if(iter->patient_ptr->disease_state==COLONIZED||iter->patient_ptr->disease_state==INFECTED){

                       if(iter->patient_ptr->disease_state==COLONIZED) {//awareness state if believed (correctly) to be positive
                           CorIstate1=KNOWN_C;
                           CorIstate2=BELIEVED_S_ACTUALLY_C;
                       } else {
                           CorIstate1=KNOWN_I;
                           CorIstate2=BELIEVED_S_ACTUALLY_I;
                       }
                       if(iter->patient_ptr->everpositiveswabthisadmission){
                        if(iter->patient_ptr->num_consec_neg_screens_following_a_pos>=policy.numbernegativeuntilconsideredcleared){
                          iter->patient_ptr->awareness_state= CorIstate2; //i.e. believed S actually C or I
                        }else {
                            iter->patient_ptr->awareness_state= CorIstate1; //Known C or I since not enough consec -ve swabs to be considered negative yet
                        }
                       } else { //patient has never been positive this admission
                        iter->patient_ptr->awareness_state= CorIstate2; //believed S actually C or I
                       }
                    } else { //patient is susceptible
                        if(iter->patient_ptr->everpositiveswabthisadmission){
                            if(iter->patient_ptr->num_consec_neg_screens_following_a_pos>=policy.numbernegativeuntilconsideredcleared){
                                iter->patient_ptr->awareness_state= KNOWN_S;
                            } else {
                                iter->patient_ptr->awareness_state=BELIEVED_C_ACTUALLY_S; //believed C since previos postive and not enough consecutive negatives
                            }
                        } else { //i.e. not positive so far on this admission
                          iter->patient_ptr->awareness_state= KNOWN_S;
                        }
                    }

                    if(iter->patient_ptr->awareness_state==KNOWN_S||iter->patient_ptr->awareness_state==BELIEVED_S_ACTUALLY_C||BELIEVED_S_ACTUALLY_I){

                        //cout<< "awareness state is   "<<int(iter->patient_ptr->awareness_state)<<"\n";

                        if(iter->patient_ptr->awareness_state!=old_awareness_state){ //so if they are believed susceptible and awareness state has changed...
                         remove_control_measures_for_those_believed_negative(iter->patient_ptr, rng); //remove any control measures in place
                        }
                    }
                }// end of iter->patient_ptr->ICU_no!=99
                iter->patient_ptr->previous_screen_result=NEGATIVE;

    //            if (iter->patient_ptr->ICU_no==99){
    //            cout<<"TRYING TO ACT ON SCREEN OF SOMEONE WHO HAS ALREADY BEEN DISCHARGED"<<"\n";
    //            // note that this has been commented out as it is possible that this happens - for example
    //            //the patient may die are be discharged in the interval between a swab being taken and the results received.
    //            }



            break;

            case SCHEDULED_SCREEN: //a planned screen for a particular patient (not part of a planned whole-ward screen)

                if  (iter->patient_ptr->ICU_no!=99 ){

                    scheduled_screening(iter->patient_ptr, rng); //calls scheduled_screening routine when time for a scheduled screen comes round

                }


            break;


            case CLINICAL_SCREEN: //occurs when a patient becomes clinically infected


                if  (iter->patient_ptr->ICU_no!=99 ){

                    clinical_screening(iter->patient_ptr, rng); //calls clinical_screening routine (after person turns I and after delaybeforeclinicalswab has passed)

                }


            break;


            case FINISH_TREATMENT:    //decolonizatino treatment ceases
            if  (iter->patient_ptr->ICU_no!=99){

                if (int(iter->event_type) == FINISH_TREATMENT){
                    end_decolonization(iter->patient_ptr, rng);
                }


            }

            } //end of switch


        }//end of for loop

      future_intervention_events.erase (time);//can delete events scheduled for current time as they have all now been dealth with



        /** ------------ DAILY OUTPUTS---------------------------------------
        */

      //      cout<<"output level is "<< OUTPUTLEVEL;
      if(OUTPUTLEVEL>0){

        cout<<"Sim "<<numsim<<" Day "<<time<<" StoC "<<StoC<<" StoI "<<StoI<<" CtoI "<<CtoI<<" acquiredInf "<<StoI+CtoI<<" discharges  " <<discharges<<" admissions "<<admissions<<" deaths "<<deaths<<" \n";
        cout<<"Sim "<<numsim<<" Day "<<time<<" app iso days "<<appisodays<< " inapp iso days "<<inappisodays<< "\n";

      }

    }//end of perform events





void patients::process_movements(gsl_rng *rng){
//see if patients die or are discharged each day - and if so schedule these events

    //cout<<"process_movements"<<"\n";

    //each day, go though all patients to see if they die or are discharged
  int num_dis_probs_sus=dis_prob_SUS.size();  //holds number of days with defined dishcarge probs for susceptibles
  //  cout<<"\n\n******dis_prob_SUS*********************************** "<< num_dis_probs_sus<<"  *********************************\n";
  int num_death_probs_sus=death_prob_SUS.size();  //holds number of days with defined death probs for susceptibles
  //cout<<"\n \n ******death prob SUS*********************************** "<< num_death_probs_sus<<"  *********************************\n";



  int num_dis_probs_inf=dis_prob_INF.size();  //holds number of days with defined dishcarge probs for infecteds
  int num_death_probs_inf=death_prob_INF.size();  //holds number of days with defined death probs for infecteds

 // int num_home_duration_sus=home_duration_SUS.size();

  int daysinhosp; //number of days in hosp so far

        for(int i=0 ; i<DEFAULT_ICU_SIZE; ++i){  //currently assume 100% bed occupancy. If lower just loop over occupied beds

            struct events nextdischargeevent;
            struct events nextdeathevent;
            struct events nextreadmissionevent;

            nextdischargeevent.event_type=ICU_DIS_AND_ADMIT;
            nextreadmissionevent.event_type=READMISSION;
            nextdeathevent.event_type=DEATH;


     //                 if (ICU_patients[i].ICU_patient_ptr->disease_state==SUSCEPTIBLE|ICU_patients[i].ICU_patient_ptr->disease_state==COLONIZED|ICU_patients[i].ICU_patient_ptr->disease_state==RECOVEREDFROMCOL|ICU_patients[i].ICU_patient_ptr->disease_state==RECOVEREDFROMINF){

            if (ICU_patients[i].ICU_patient_ptr->disease_state==SUSCEPTIBLE||ICU_patients[i].ICU_patient_ptr->disease_state==COLONIZED){

                //create arrays of discharge and death probabilities for all SUS / COL people in ICU on each day
                double dischargeSUS [DEFAULT_ICU_SIZE];
                double deathSUS [DEFAULT_ICU_SIZE];

       daysinhosp=ICU_patients[i].ICU_patient_ptr->no_days_in_hos;
       if(daysinhosp >= num_dis_probs_sus){ //if day of stay is more than number of days we have dis probs for...
         dischargeSUS[i]=dis_prob_SUS[num_dis_probs_sus-1];
       } else {
         dischargeSUS[i] = dis_prob_SUS[ICU_patients[i].ICU_patient_ptr->no_days_in_hos];
       }

       if(daysinhosp >= num_death_probs_sus){
         deathSUS[i]=death_prob_SUS[num_death_probs_sus-1];
       } else {
         deathSUS[i] = death_prob_SUS[ICU_patients[i].ICU_patient_ptr->no_days_in_hos];
       }


                float randnum;
                randnum = gsl_rng_uniform(rng);


                //cout<<dischargeSUS[i]<<"discharge prob for person" <<i<<"\n";

                //does discharge occur?
                if (randnum < dischargeSUS[i]){

                                // cout<< "DISCAHRGE Q'd for pateint "<<ICU_patients[i].ICU_patient_ptr->ICU_no<<"\n";

                                //queue up discharge
                    nextdischargeevent.patient_ptr = ICU_patients[i].ICU_patient_ptr;
                    future_movement_events[time].push_back(nextdischargeevent);

               // if (randnum < readmissionSUS[i]){

                 //   nextreadmissionevent.patient_ptr = ICU_patients[i].ICU_patient_ptr;
                   // future_movement_events[time].push_back(nextreadmissionevent);

//                }


                } else {
                    if (randnum-dischargeSUS[i]<deathSUS[i]){


                        // cout<< "DEATH Q'd for pateint "<<ICU_patients[i].ICU_patient_ptr->ICU_no<<"\n";

                        //queue up death
                        nextdeathevent.patient_ptr = ICU_patients[i].ICU_patient_ptr;
                        future_movement_events[time].push_back(nextdeathevent);
                    }
                }
            }//if sus or col


    //        if (ICU_patients[i].ICU_patient_ptr->disease_state==INFECTED|ICU_patients[i].ICU_patient_ptr->disease_state==INFECTEDFROMCOL){
        if (ICU_patients[i].ICU_patient_ptr->disease_state==INFECTED){
                //create an array of discharge and death probabilities for all INF people in ICU on each day
                double dischargeINF [DEFAULT_ICU_SIZE];//this is too big though - doesnt need to be whole size of isolation ward
                double deathINF [DEFAULT_ICU_SIZE];

       daysinhosp=ICU_patients[i].ICU_patient_ptr->no_days_in_hos;
     //  if(daysinhosp >= num_dis_probs_inf){

         if(daysinhosp >= num_dis_probs_inf){
         dischargeINF[i]=dis_prob_INF[num_dis_probs_inf-1];
       } else {
         dischargeINF[i] = dis_prob_INF[ICU_patients[i].ICU_patient_ptr->no_days_in_hos];
       }

       if(daysinhosp >= num_death_probs_inf){
         deathINF[i]=death_prob_INF[num_death_probs_inf-1];
       } else {
         deathINF[i] = death_prob_INF[ICU_patients[i].ICU_patient_ptr->no_days_in_hos];
       }


                float randnum;
                randnum = gsl_rng_uniform(rng);

                //cout<<dischargeINF[i]<<"discharge prob for person" <<i<<"\n";

                //does discharge occur?
                if (randnum < dischargeINF[i]){

                // cout<< "DISCAHRGE Q'd for pateint "<<ICU_patients[i].ICU_patient_ptr->ICU_no<<"\n";

                //queue up discharge
                    nextdischargeevent.patient_ptr = ICU_patients[i].ICU_patient_ptr;
                    future_movement_events[time].push_back(nextdischargeevent);



                } else {
                    if (randnum-dischargeINF[i]<deathINF[i]){


                        //cout<< "DEATH Q'd for pateint "<<ICU_patients[i].ICU_patient_ptr->ICU_no<<"\n";

                        //queue up death
                        nextdeathevent.patient_ptr = ICU_patients[i].ICU_patient_ptr;
                        future_movement_events[time].push_back(nextdeathevent);
                    }
                }


            }//if infected

        }//for all patients

}//process_movements







void patients::process_infections(gsl_rng *rng){

  //  ofstream prev_file("infections_over_time.txt",  ios::app );



  if (DEFAULT_DEBUG==32){

  cout<<"*****PROCESS INFECTIONS*****"<<"\n";
  }

  //function to:
  //count up numbers of S C and I
  //loop through all sus and col and calculate their probabilities of becoming colonized/infected
  //randomly choose if these events occur based on the calcualted probabilies
  //if they do, add them to the event queue

  double probcolgiven1inf ;
  double probinfgiven1inf;
  double probcol;
  double probinf;
  double progprob;

        int numberinfected=0;
        int numbercolonized=0;
        int numbersusceptible=0;
        int numbercolorinf = 0;
        double C_IQ_count = 0;  //total colonization pressure due to colonized patients
        double I_IQ_count = 0;  //total colonization pressure due to infected patients
        double total_IQ_in_ward = 0; //total colonization pressure due to infecteds and those colonized
 float randnum1, randnum2;

  for(int i=0 ; i<DEFAULT_ICU_SIZE; ++i){
    //   if(ICU_patients[i].ICU_patient_ptr->disease_state==SUSCEPTIBLE|ICU_patients[i].ICU_patient_ptr->disease_state==RECOVEREDFROMCOL|ICU_patients[i].ICU_patient_ptr->disease_state==RECOVEREDFROMINF){ - Ben - not sure this makes sense...disease state should be susceptilbe, infected or colonized

    if(ICU_patients[i].ICU_patient_ptr->disease_state==SUSCEPTIBLE){
    ++numbersusceptible;
   }


   if(ICU_patients[i].ICU_patient_ptr->disease_state==COLONIZED){
    //colonized[numbercolonized]=i;

    ++numbercolonized;
    C_IQ_count = C_IQ_count + double (ICU_patients[i].ICU_patient_ptr->IQ);

    //cout<<ICU_patients[i].ICU_patient_ptr->patientid<<"  colonized patient ID"<<"\n";


   }

  // prev_file<<numbersusceptible<<" \t"<<numbercolonized<<" \t"<<numberinfected<<" \n";


   //         if(ICU_patients[i].ICU_patient_ptr->disease_state==INFECTED|ICU_patients[i].ICU_patient_ptr->disease_state==INFECTEDFROMCOL ){  - Ben - dont' see why we need an INFECTEDFROMCOL disease state
         if(ICU_patients[i].ICU_patient_ptr->disease_state==INFECTED){
    //infecteds[numberinfected]=i;
    ++numberinfected;
    I_IQ_count = I_IQ_count + double (ICU_patients[i].ICU_patient_ptr->IQ);

    //cout<<ICU_patients[i].ICU_patient_ptr->patientid<<"  infected patient ID"<<"\n";
    //cout<<double (ICU_patients[i].ICU_patient_ptr->IQ)<<"IQ of patient"<<"\n";
    //cout<<I_IQ_count<<"  I_IQ_count ********************************"<<"\n";

            }

                //cout<<int(ICU_patients[i].ICU_patient_ptr->disease_state)<< "disease states of ICU population"<<"\n";
                //cout<<int(ICU_patients[i].ICU_patient_ptr->ICU_no)<< "ICU_no"<<"\n";



    //          if(ICU_patients[i].ICU_patient_ptr->disease_state==COLANDINF  ){
    //    //infecteds[numberinfected]=i;
    //    ++numbercolandinf;
    //    //cout<<numbercolandinf<< " number infected and colonized" << "\n";
    //   }

    }//end of for (looping through all patients in the ICU)

            /**
    if (DEFAULT_DEBUG==32){
                    cout<<numbersusceptible<< " number susceptible" << "\n";
    }


       if (DEFAULT_DEBUG==32){
                    cout<<numbercolonized<< " number colonized" << "\n";
    }


                if (DEFAULT_DEBUG==32){
                    cout<<numberinfected<< " number infected" << "\n";
    }
            */
        numbercolorinf = numbercolonized+numberinfected ;

        //cout<<numbercolorinf <<" = number col or inf" << "\n";
        total_IQ_in_ward = C_IQ_count+I_IQ_count;
 //        cout<<total_IQ_in_ward<< "total IQ in ward"<<"\n";

        if (DEFAULT_DEBUG==32){
          cout<<numbercolorinf <<" = number col or inf" << "\n";
        }



////////////////////////////////    CALCULATION OF ALL INF AND COL EVENTS  /////////////////////

    struct events nextcolevent = {COLONIZE, NULL};
    struct events nextinfevent = {INFECT, NULL};
    struct events nextinffromcolevent = {INFECT, NULL};
    //struct events nextcolandinfevent = {COLANDINF, NULL};


            //loop through ICU and calcualte prob of all individuals becoming col or inf
            for (int m=0; m<DEFAULT_ICU_SIZE; m++){

                    //  do any sus become col or inf?

       // again - don't see that we need all these disaease states                    if (ICU_patients[m].ICU_patient_ptr->disease_state==SUSCEPTIBLE|ICU_patients[m].ICU_patient_ptr->disease_state==RECOVEREDFROMCOL|ICU_patients[m].ICU_patient_ptr->disease_state==RECOVEREDFROMINF){

                  if (ICU_patients[m].ICU_patient_ptr->disease_state==SUSCEPTIBLE){


      //section below commented out - no need for an array which is repeatedly declared.
//           double probcolgiven1inf [DEFAULT_ICU_SIZE];
//                         double probinfgiven1inf [DEFAULT_ICU_SIZE];
//                         double probcol [DEFAULT_ICU_SIZE];
//                         double probinf [DEFAULT_ICU_SIZE];

//                         probcolgiven1inf[m] = ICU_patients[m].ICU_patient_ptr->Pdc; //find susceptibilty of patient m to col
//                         probinfgiven1inf[m] = ICU_patients[m].ICU_patient_ptr->Pdi; //find susceptibilty of patient m to inf
//                         probcol[m] = 1-(pow((1-probcolgiven1inf[m]), total_IQ_in_ward));//calculate prob of patient m becoming col (given there are number colorinf infectious individuals  - assuming col and inf transmit equally)
//                         probinf[m] = 1-(pow((1-probinfgiven1inf[m]), total_IQ_in_ward)); //calculate prob of patient m becoming inf
//                         //cout<<total_IQ_in_ward<<"total IQ in ward" <<"\n";

                        probcolgiven1inf = ICU_patients[m].ICU_patient_ptr->Pdc; //find susceptibilty of patient m to col
                        probinfgiven1inf = ICU_patients[m].ICU_patient_ptr->Pdi; //find susceptibilty of patient m to inf
                        probcol = 1-(pow((1-probcolgiven1inf), total_IQ_in_ward));//calculate prob of patient m becoming col (given there are number colorinf infectious individuals  - assuming col and inf transmit equally)
                        probinf = 1-(pow((1-probinfgiven1inf), total_IQ_in_ward)); //calculate prob of patient m becoming inf


                        if (DEFAULT_DEBUG==32){
                            cout<<probcol<< "= prob col      " <<probinf<<" = prob inf"   << "\n";
                        }

                        randnum1 = gsl_rng_uniform(rng);

                        //does col occur?
                        if (randnum1 < probcol ){//probinf[m]<

                            nextcolevent.patient_ptr=ICU_patients[m].ICU_patient_ptr;
                            future_infection_events[time+1].push_back(nextcolevent);

                            if (DEFAULT_DEBUG==32){
                                cout<< " queued up a colonization to person  "<<nextcolevent.patient_ptr->patientid<<"\n";
                            }

                        }

                        randnum1 = gsl_rng_uniform(rng);
                        //does inf oocur
                        if (randnum1 < probinf){

                            nextinfevent.patient_ptr=ICU_patients[m].ICU_patient_ptr;
                            future_infection_events[time+1].push_back(nextinfevent);


                            if (DEFAULT_DEBUG==32){
                                cout<< " queued up an infection to pateint    "<<nextinfevent.patient_ptr->patientid<< "   who was  "<<int(nextinfevent.patient_ptr->disease_state)<< "\n";
                            }
                        }

                    }//if susceptible

                    //do any colonized progress to inf?

                    if (ICU_patients[m].ICU_patient_ptr->disease_state==COLONIZED){



                        randnum1 = gsl_rng_uniform(rng);
                        //does inf oocur
                        if (randnum1 < probinf){

                            nextinffromcolevent.patient_ptr=ICU_patients[m].ICU_patient_ptr;
                            future_infection_events[time+1].push_back(nextinffromcolevent);


                            if (DEFAULT_DEBUG==32){
                              //  cout<< " queued up an inf from col event (transmission)    "<<nextinffromcolevent.patient_ptr->patientid<< "   who was  "<<int(nextinffromcolevent.patient_ptr->disease_state)<< "\n";
                            }
                        }


                        progprob = ICU_patients[m].ICU_patient_ptr->ProgProb;


                        randnum2 = gsl_rng_uniform(rng);

                        if (randnum2 < progprob){

                          nextinffromcolevent.patient_ptr=ICU_patients[m].ICU_patient_ptr;
                          future_infection_events[time+1].push_back(nextinffromcolevent);


                          if (DEFAULT_DEBUG==32){
                            cout<< " sceduled inf from col event (progression) to patient    "<<nextinffromcolevent.patient_ptr->patientid<<"\n";
                          }

                        }

                    }//if colonized

            }//for (loop through ICU)

}//end of process_infections


void patients::app_or_inapp_isolation_days(){

  for(int i=0 ; i<DEFAULT_ICU_SIZE; ++i){

      if(ICU_patients[i].ICU_patient_ptr!=NULL && ICU_patients[i].ICU_patient_ptr->isolation==true){

      //cout<<"patient is in isoaltion and their disease state is "<<int(ICU_patients[i].ICU_patient_ptr->disease_state)<<"\n";

                //      Appropriate isolation days-----------------------------------------------

          if(ICU_patients[i].ICU_patient_ptr->disease_state==COLONIZED || ICU_patients[i].ICU_patient_ptr->disease_state==INFECTED ){
                    ++appisodays;
                    ++cumulative_appisodays;
          }

           //      Inappropriate isolation days --------------------------------------------

                //      [put in isolation but shouldnt be]

                if(ICU_patients[i].ICU_patient_ptr->disease_state==SUSCEPTIBLE ){
                    ++inappisodays;
                    ++cumulative_inappisodays;


          }
      }

            if(ICU_patients[i].ICU_patient_ptr!=NULL && ICU_patients[i].ICU_patient_ptr->isolation==false){

                    if(ICU_patients[i].ICU_patient_ptr->disease_state==COLONIZED || ICU_patients[i].ICU_patient_ptr->disease_state==INFECTED ){
                    ++unisodays;
                    ++cumulative_unisodays;


//cout<<int(ICU_patients[i].ICU_patient_ptr->disease_state)<<"   Disease state of unisolated patient " <<"\n";
//cout<<ICU_patients[i].ICU_patient_ptr->patientid <<" id number of inisolated patient "<<"\n";


                    }

            }


//cout<<appisodays <<" appisodays "<<"\n";
//cout<<cumulative_appisodays <<" cumulativive_appisodays "<<"\n";


//cout<<inappisodays <<" inappisodays "<<"\n";
//cout<<cumulative_inappisodays <<" cumulativive_inappisodays "<<"\n";

//cout<<unisodays <<" unisodays "<<"\n";
//cout<<cumulative_unisodays <<" cumulativive_unisodays "<<"\n";

  }

}






void patients::bed_day_counts(){

  for(int i=0 ; i<DEFAULT_ICU_SIZE; ++i){
 //  if(ICU_patients[i].ICU_patient_ptr->disease_state==SUSCEPTIBLE|ICU_patients[i].ICU_patient_ptr->disease_state==RECOVEREDFROMCOL|ICU_patients[i].ICU_patient_ptr->disease_state==RECOVEREDFROMINF){
  if(ICU_patients[i].ICU_patient_ptr->disease_state==SUSCEPTIBLE){
    ++susbedday;
   }
   if(ICU_patients[i].ICU_patient_ptr->disease_state==COLONIZED){
    ++colbedday;
   }
//   if(ICU_patients[i].ICU_patient_ptr->disease_state==INFECTED|ICU_patients[i].ICU_patient_ptr->disease_state==INFECTEDFROMCOL){
   if(ICU_patients[i].ICU_patient_ptr->disease_state==INFECTED){
    ++infbedday;
   }

   if(ICU_patients[i].ICU_patient_ptr->isolation==true){
    ++isobedday;
   }


  }




} //end of bed_day_counts

//each day while someone is in the ICU adds one to their no_days_in_hos
void patients::LOS_counts(){

  for(int i=0 ; i<DEFAULT_ICU_SIZE; ++i){
                    ICU_patients[i].ICU_patient_ptr->no_days_in_hos=ICU_patients[i].ICU_patient_ptr->no_days_in_hos + 1;
  }
} //end of LOS_counts




void patients::setpolicy(gsl_rng *rng, const int policynumber){
 //sets the current screening and intervention policy and screening and isolation parameters
 //this routine will be modified so that policy can be read from the command line, env variables or interactively from the user
 //but at the moment it just chooses from predefined policies according to the policynumber
  // cout<<"setpolicy not currently  implemented";


  // sentivity and spcificity of different cultures are sampled randomly based on global means and SDs (which can be set on the commmand line)

  int isolation_capacity=get_isocap();
  //cout<<"isolation_capacity="<<isolation_capacity<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
  double primary_iso_effectiveness =  get_effect_of_ISO(rng);
  double secondary_iso_effectiveness=  get_effect_of_secISO(rng);

  double bwash_pdc_effectiveness = get_pdc_effect_of_bodywash(rng);
  double bwash_pdi_effectiveness = get_pdi_effect_of_bodywash(rng);
//  double bwash_IQ_effectiveness = get_IQ_effect_of_bodywash(rng);
  double bwash_progprob_effectiveness = get_progprob_effect_of_bodywash(rng);


 conventional_culture.sensitivity= get_CC_sensitivity(rng);;
 conventional_culture.specificity= get_CC_specificity(rng);
 conventional_culture.tat=4;

 chromagar.sensitivity=get_CA_sensitivity(rng);
 chromagar.specificity=get_CA_specificity(rng);
 chromagar.tat=3;

 chromagar_early.sensitivity=get_CA_early_sensitivity(rng);
 chromagar_early.specificity=get_CA_early_specificity(rng);
 chromagar_early.tat=2;// 16-24h for test + extra bits

 pcr.sensitivity= get_PCR_sensitivity(rng);//0.88; //for IDI-MRSA
 pcr.specificity=get_PCR_specificity(rng);//0.84;//for IDI-MRSA
 pcr.tat=1;

 ideal.sensitivity= 0.884; //characteristics of an ideal test (with best propoerties of all available technologies).
 ideal.specificity=0.9713;
 ideal.tat=1;

  struct screeningpolicycomponent convculture1; //conventional culture non-targeted
  convculture1.technology=CC; //other options are CA (chromagar) CA_early and PCR
  convculture1.targeted=false;
  convculture1.previouslypositive=false;

  struct screeningpolicycomponent ca1; //CA non-targeted
  ca1.technology=CA;
  ca1.targeted=false;
  ca1.previouslypositive=false;

  struct screeningpolicycomponent ca_early1; //CA non-targeted
  ca_early1.technology=CA_early;
  ca_early1.targeted=false;
  ca_early1.previouslypositive=false;


  struct screeningpolicycomponent pcr1; //PCR non-targeted
  pcr1.technology=PCR;
  pcr1.targeted=false;
  pcr1.previouslypositive=false;

  struct screeningpolicycomponent ideal1;
  ideal1.technology=I;  //I for ideal
  ideal1.targeted=false;
  ideal1.previouslypositive=false;

 struct  interventionpolicy  policy1; //no screening, just transmission plus clinical cultures from infecteds and contact precuations for those

  policy1.screeninghighriskonly = false;

  policy1.clinicalscreening.push_back(convculture1);
  policy1.preemptiveisolationforall=false;
  policy1.preemptiveisolationforhighrisk=false;
  policy1.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy1.isolateifpositive=true;
  policy1.blanketdecolforall=false;
  policy1.decolforhighrisk=false;
  policy1.decolifpositive=false;
  policy1.timesteptoimplement=600;//i.e. at what timestep does the policy kick in
  policy1.proportionscreenedonadmission=0;
  policy1.proportionscreenedweekly=0;
  policy1.proportionscreenedondischarge=0;
  policy1.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken


    //  ===================================JUST TRANS======================

  //  ===================================admission screening and weekly screening on a monday with conv culture======================
    //screening with conv culture accompanied by isolation of positives


  struct  interventionpolicy  policy2; //discharge and admission screening and weekly screening on a monday


  policy2.screeninghighriskonly = false;
  policy2.admissionscreening.push_back(convculture1);
  policy2.dischargescreening.push_back(convculture1);
  policy2.weekdayscreening[0].push_back(convculture1);
  policy2.clinicalscreening.push_back(convculture1);
  policy2.preemptiveisolationforall=false;
  policy2.preemptiveisolationforhighrisk=false;
  policy2.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy2.isolateifpositive=true;
  policy2.blanketdecolforall=false;
  policy2.decolforhighrisk=false;
  policy2.decolifpositive=false;
  policy2.timesteptoimplement=365; //365;//i.e. at what timestep does the policy kick in
  policy2.proportionscreenedonadmission=1;
  policy2.proportionscreenedweekly=1;
  policy2.proportionscreenedondischarge=0;

  policy2.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken


  //  ===================================admission screening and weekly screening on a monday with  chromogenic agar======================
 //   contact isolation of positives



  struct  interventionpolicy  policy3; //admission screening and weekly screening on a monday with ca====================

  policy3.screeninghighriskonly = false;
  policy3.admissionscreening.push_back(ca1);
  policy3.weekdayscreening[0].push_back(ca1);
  policy3.clinicalscreening.push_back(ca1);
// policy3.dischargescreening.push_back(convculture1);
  //currently no weeklypostadmission or discharge screening
  policy3.preemptiveisolationforall=false;
  policy3.preemptiveisolationforhighrisk=false;
  policy3.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy3.isolateifpositive=true;
  policy3.blanketdecolforall=false;
  policy3.decolforhighrisk=false;
  policy3.decolifpositive=false;
  policy3.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy3.proportionscreenedonadmission=1;
  policy3.proportionscreenedweekly=1;
  policy3.proportionscreenedondischarge=0;
  policy3.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken

//  ===================================admission screening and weekly screening on a monday with  chromogenic agar (taking early result)======================
   //+ contact isolation of positives

  struct  interventionpolicy  policy4; //admission screening and weekly screening on a monday with ca_early====================
                                        // and amended with the ca 48 hour result

  policy4.screeninghighriskonly = false;
  policy4.admissionscreening.push_back(ca_early1);
  policy4.weekdayscreening[0].push_back(ca_early1);
  policy4.clinicalscreening.push_back(ca_early1);
  policy4.admissionscreening.push_back(ca1);
  policy4.weekdayscreening[0].push_back(ca1);
  policy4.clinicalscreening.push_back(ca1);

  //currently no weeklypostadmission or discharge screening
  policy4.preemptiveisolationforall=false;
  policy4.preemptiveisolationforhighrisk=false;
  policy4.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy4.isolateifpositive=true;
  policy4.blanketdecolforall=false;
  policy4.decolforhighrisk=false;
  policy4.decolifpositive=false;
  policy4.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy4.proportionscreenedonadmission=1;
  policy4.proportionscreenedweekly=1;
  policy4.proportionscreenedondischarge=0;
  policy4.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken

 //  ===================================admission screening and weekly screening on a monday with  PCR======================
    //screening with PCR accompanied by isolation of positives


  struct  interventionpolicy  policy5; //admission screening and weekly screening on a monday with pcr====================

  policy5.screeninghighriskonly = false;
  policy5.admissionscreening.push_back(pcr1);
  policy5.weekdayscreening[0].push_back(pcr1);
  policy5.clinicalscreening.push_back(pcr1);
  //currently no weeklypostadmission or discharge screening
  policy5.preemptiveisolationforall=false;
  policy5.preemptiveisolationforhighrisk=false;
  policy5.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy5.isolateifpositive=true;
  policy5.blanketdecolforall=false;
  policy5.decolforhighrisk=false;
  policy5.decolifpositive=false;
  policy5.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy5.proportionscreenedonadmission=1;
  policy5.proportionscreenedweekly=1;
  policy5.proportionscreenedondischarge=0;
  policy5.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken


//  ===================================admission screening and weekly screening on a monday with  CC======================
//===================================== + pre-emptive isolation for high risk and isolation of known positives when screening results come in===========================================================

  struct  interventionpolicy  policy6;

  policy6.screeninghighriskonly = false;
  policy6.admissionscreening.push_back(convculture1);
  policy6.weekdayscreening[0].push_back(convculture1);
  policy6.clinicalscreening.push_back(convculture1);
  //currently no weeklypostadmission or discharge screening
  policy6.preemptiveisolationforall=false;
  policy6.preemptiveisolationforhighrisk=true;
  policy6.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy6.isolateifpositive=true;
  policy6.blanketdecolforall=false;
  policy6.decolforhighrisk=false;
  policy6.decolifpositive=false;
  policy6.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy6.proportionscreenedonadmission=1;
  policy6.proportionscreenedweekly=1;
  policy6.proportionscreenedondischarge=0;
  policy6.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken


//  ===================================no screening======================
//===================================== just blanket isolation for all ===========================================================

  struct  interventionpolicy  policy7;


  policy7.screeninghighriskonly = false;
  policy7.preemptiveisolationforall=true;
  policy7.preemptiveisolationforhighrisk=false;
  policy7.clinicalscreening.push_back(convculture1);
  policy7.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy7.isolateifpositive=false;
  policy7.blanketdecolforall=false;
  policy7.decolforhighrisk=false;
  policy7.decolifpositive=false;
  policy7.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy7.proportionscreenedonadmission=1;
  policy7.proportionscreenedweekly=1;
  policy7.proportionscreenedondischarge=0;
  policy7.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken


//  ===================================no screening======================
//===================================== just blanket isolation for high risk ===========================================================

//(with this policy also need contact precuations from known positives, highlighted through clinical cultures)

  struct  interventionpolicy  policy8;

  policy8.screeninghighriskonly = false;
  policy8.clinicalscreening.push_back(convculture1);

  policy8.preemptiveisolationforall=false;
  policy8.preemptiveisolationforhighrisk=true;
  policy8.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy8.isolateifpositive=true;
  policy8.blanketdecolforall=false;
  policy8.decolforhighrisk=false;
  policy8.decolifpositive=false;
  policy8.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy8.proportionscreenedonadmission=1;
  policy8.proportionscreenedweekly=1;
  policy8.proportionscreenedondischarge=0;
  policy8.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken



//  ===================================screening high risk only (with  CC)======================
//===================================== + contact precautions for those who are positive ===========================================================

 struct  interventionpolicy  policy9;

  policy9.screeninghighriskonly = true;
  policy9.admissionscreening.push_back(convculture1);

  //NOT SURE THIS IS WORKING

  policy9.weekdayscreening[0].push_back(convculture1);
  policy9.clinicalscreening.push_back(convculture1);

  //currently no weeklypostadmission or discharge screening
  policy9.preemptiveisolationforall=false;
  policy9.preemptiveisolationforhighrisk=false;
  policy9.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy9.isolateifpositive=true;
  policy9.blanketdecolforall=false;
  policy9.decolforhighrisk=false;
  policy9.decolifpositive=false;
  policy9.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy9.proportionscreenedonadmission=1;
  policy9.proportionscreenedweekly=1;
  policy9.proportionscreenedondischarge=0;
  policy9.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken



//  ===================================screening high risk only (with  CA)======================
//===================================== + contact precautions for those who are positive ===========================================================

  struct  interventionpolicy  policy10;

  policy10.screeninghighriskonly = true;
  policy10.admissionscreening.push_back(ca1);

  //NOT SURE THIS IS WORKING

  policy10.weekdayscreening[0].push_back(ca1);
  policy10.clinicalscreening.push_back(ca1);

  //currently no weeklypostadmission or discharge screening
  policy10.preemptiveisolationforall=false;
  policy10.preemptiveisolationforhighrisk=false;
  policy10.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy10.isolateifpositive=true;
  policy10.blanketdecolforall=false;
  policy10.decolforhighrisk=false;
  policy10.decolifpositive=false;
  policy10.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy10.proportionscreenedonadmission=1;
  policy10.proportionscreenedweekly=1;
  policy10.proportionscreenedondischarge=0;
  policy10.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken


//  ===================================screening high risk only (with  CA_early)======================
//===================================== + contact precautions for those who are positive ===========================================================

  struct  interventionpolicy  policy11;

  policy11.screeninghighriskonly = true;
  policy11.admissionscreening.push_back(ca_early1);

  //NOT SURE THIS IS WORKING

  policy11.weekdayscreening[0].push_back(ca_early1);
  policy11.clinicalscreening.push_back(ca_early1);

  //currently no weeklypostadmission or discharge screening
  policy11.preemptiveisolationforall=false;
  policy11.preemptiveisolationforhighrisk=false;
  policy11.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy11.isolateifpositive=true;
  policy11.blanketdecolforall=false;
  policy11.decolforhighrisk=false;
  policy11.decolifpositive=false;
  policy11.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy11.proportionscreenedonadmission=1;
  policy11.proportionscreenedweekly=1;
  policy11.proportionscreenedondischarge=0;
  policy11.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken


//  ===================================screening high risk only (with  PCR)======================
//===================================== + contact precautions for those who are positive ===========================================================

  struct  interventionpolicy  policy12;

  policy12.screeninghighriskonly = true;
  policy12.admissionscreening.push_back(pcr1);

  //NOT SURE THIS IS WORKING

  policy12.weekdayscreening[0].push_back(pcr1);
  policy12.clinicalscreening.push_back(pcr1);

  //currently no weeklypostadmission or discharge screening
  policy12.preemptiveisolationforall=false;
  policy12.preemptiveisolationforhighrisk=false;
  policy12.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy12.isolateifpositive=true;
  policy12.blanketdecolforall=false;
  policy12.decolforhighrisk=false;
  policy12.decolifpositive=false;
  policy12.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy12.proportionscreenedonadmission=1;
  policy12.proportionscreenedweekly=1;
  policy12.proportionscreenedondischarge=0;
  policy12.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken


 //  ===================================admission screening and weekly screening on a monday with  an ideal test======================
    //screening with an ideal test accompanied by isolation of positives

  struct  interventionpolicy  policy13; //admission screening and weekly screening on a monday with pcr====================


  policy13.screeninghighriskonly = false;
  policy13.admissionscreening.push_back(ideal1);
  policy13.weekdayscreening[0].push_back(ideal1);
  policy13.clinicalscreening.push_back(ideal1);
  //currently no weeklypostadmission or discharge screening
  policy13.preemptiveisolationforall=false;
  policy13.preemptiveisolationforhighrisk=false;
  policy13.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy13.isolateifpositive=true;
  policy13.blanketdecolforall=false;
  policy13.decolforhighrisk=false;
  policy13.decolifpositive=false;
  policy13.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy13.proportionscreenedonadmission=1;
  policy13.proportionscreenedweekly=1;
  policy13.proportionscreenedondischarge=0;
  policy13.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken


//  ===================================screening high risk only (with  ideal test)======================
//===================================== + contact precautions for those who are positive ===========================================================

  struct  interventionpolicy  policy14;

  policy14.screeninghighriskonly = true;
  policy14.admissionscreening.push_back(ideal1);

  //NOT SURE THIS IS WORKING

  policy14.weekdayscreening[0].push_back(ideal1);
  policy14.clinicalscreening.push_back(ideal1);

  //currently no weeklypostadmission or discharge screening
  policy14.preemptiveisolationforall=false;
  policy14.preemptiveisolationforhighrisk=false;
  policy14.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy14.isolateifpositive=true;
  policy14.blanketdecolforall=false;
  policy14.decolforhighrisk=false;
  policy14.decolifpositive=false;
  policy14.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy14.proportionscreenedonadmission=1;
  policy14.proportionscreenedweekly=1;
  policy14.proportionscreenedondischarge=0;
  policy14.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken



//  ===================================no screening=====================
//===================================== + blanket decolonisation for all ===========================================================

  struct  interventionpolicy  policy15;

  policy15.screeninghighriskonly = false;
  policy15.clinicalscreening.push_back(convculture1);

  policy15.preemptiveisolationforall=false;
  policy15.preemptiveisolationforhighrisk=false;
  policy15.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy15.isolateifpositive=false;
  policy15.blanketdecolforall=true;
  policy15.decolforhighrisk=false;
  policy15.decolifpositive=false;
  policy15.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy15.proportionscreenedonadmission=1;
  policy15.proportionscreenedweekly=1;
  policy15.proportionscreenedondischarge=0;
  policy15.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken


//  ===================================no screening=====================
//===================================== + decolonisation for high risk ===========================================================

//(this still needs clinical swabs - with decolonisation of those found positive)

  struct  interventionpolicy  policy16;

  policy16.screeninghighriskonly = false;
  policy16.clinicalscreening.push_back(convculture1);

  policy16.preemptiveisolationforall=false;
  policy16.preemptiveisolationforhighrisk=false;
  policy16.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy16.isolateifpositive=false;
  policy16.blanketdecolforall=false;
  policy16.decolforhighrisk=true;
  policy16.decolifpositive=true;
  policy16.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy16.proportionscreenedonadmission=1;
  policy16.proportionscreenedweekly=1;
  policy16.proportionscreenedondischarge=0;
  policy16.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken



  //  ===================================admission screening and weekly screening on a monday with  conventional culture======================
//+ decolonisation of positives

  struct  interventionpolicy  policy17;

  policy17.screeninghighriskonly = false;
  policy17.admissionscreening.push_back(convculture1);
  policy17.weekdayscreening[0].push_back(convculture1);
  policy17.clinicalscreening.push_back(convculture1);
  //currently no weeklypostadmission or discharge screening
  policy17.preemptiveisolationforall=false;
  policy17.preemptiveisolationforhighrisk=false;
  policy17.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy17.isolateifpositive=false;
  policy17.blanketdecolforall=false;
  policy17.decolforhighrisk=false;
  policy17.decolifpositive=true;
  policy17.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy17.proportionscreenedonadmission=1;
  policy17.proportionscreenedweekly=1;
  policy17.proportionscreenedondischarge=0;
  policy17.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken


//  ===================================screening high risk only (with  CC)======================
//===================================== + decolonisation for those who are positive ===========================================================

  struct  interventionpolicy  policy18;

  policy18.screeninghighriskonly = true;
  policy18.admissionscreening.push_back(convculture1);

  //NOT SURE THIS IS WORKING

  policy18.weekdayscreening[0].push_back(convculture1);
  policy18.clinicalscreening.push_back(convculture1);

  //currently no weeklypostadmission or discharge screening
  policy18.preemptiveisolationforall=false;
  policy18.preemptiveisolationforhighrisk=false;
  policy18.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy18.isolateifpositive=false;
  policy18.blanketdecolforall=false;
  policy18.decolforhighrisk=false;
  policy18.decolifpositive=true;
  policy18.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy18.proportionscreenedonadmission=1;
  policy18.proportionscreenedweekly=1;
  policy18.proportionscreenedondischarge=0;
  policy18.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken



  //  ===================================admission screening and weekly screening on a monday with  chromogenic agar======================
//+ decolonisation of positives

  struct  interventionpolicy  policy19;


  policy19.screeninghighriskonly = false;
  policy19.admissionscreening.push_back(ca1);
  policy19.weekdayscreening[0].push_back(ca1);
  policy19.clinicalscreening.push_back(ca1);
  //currently no weeklypostadmission or discharge screening
  policy19.preemptiveisolationforall=false;
  policy19.preemptiveisolationforhighrisk=false;
  policy19.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy19.isolateifpositive=false;
  policy19.blanketdecolforall=false;
  policy19.decolforhighrisk=false;
  policy19.decolifpositive=true;
  policy19.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy19.proportionscreenedonadmission=1;
  policy19.proportionscreenedweekly=1;
  policy19.proportionscreenedondischarge=0;
  policy19.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken



//  ===================================screening high risk only (with  CA)======================
//===================================== + decolonisation for those who are positive ===========================================================

  struct  interventionpolicy  policy20;

  policy20.screeninghighriskonly = true;
  policy20.admissionscreening.push_back(ca1);

  //NOT SURE THIS IS WORKING

  policy20.weekdayscreening[0].push_back(ca1);
  policy20.clinicalscreening.push_back(ca1);

  //currently no weeklypostadmission or discharge screening
  policy20.preemptiveisolationforall=false;
  policy20.preemptiveisolationforhighrisk=false;
  policy20.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy20.isolateifpositive=false;
  policy20.blanketdecolforall=false;
  policy20.decolforhighrisk=false;
  policy20.decolifpositive=true;
  policy20.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy20.proportionscreenedonadmission=1;
  policy20.proportionscreenedweekly=1;
  policy20.proportionscreenedondischarge=0;
  policy20.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken




  //  ===================================admission screening and weekly screening on a monday with  PCR======================
//+ decolonisation of positives

  struct  interventionpolicy  policy21;

  policy21.screeninghighriskonly = false;
  policy21.admissionscreening.push_back(pcr1);
  policy21.weekdayscreening[0].push_back(pcr1);
  policy21.clinicalscreening.push_back(pcr1);
  //currently no weeklypostadmission or discharge screening
  policy21.preemptiveisolationforall=false;
  policy21.preemptiveisolationforhighrisk=false;
  policy21.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy21.isolateifpositive=false;
  policy21.blanketdecolforall=false;
  policy21.decolforhighrisk=false;
  policy21.decolifpositive=true;
  policy21.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy21.proportionscreenedonadmission=1;
  policy21.proportionscreenedweekly=1;
  policy21.proportionscreenedondischarge=0;
  policy21.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken



//  ===================================screening high risk only (with  PCR)======================
//===================================== + decolonisation for those who are positive ===========================================================

  struct  interventionpolicy  policy22;

  policy22.screeninghighriskonly = true;
  policy22.admissionscreening.push_back(pcr1);

  //NOT SURE THIS IS WORKING

  policy22.weekdayscreening[0].push_back(pcr1);
  policy22.clinicalscreening.push_back(pcr1);

  //currently no weeklypostadmission or discharge screening
  policy22.preemptiveisolationforall=false;
  policy22.preemptiveisolationforhighrisk=false;
  policy22.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy22.isolateifpositive=false;
  policy22.blanketdecolforall=false;
  policy22.decolforhighrisk=false;
  policy22.decolifpositive=true;
  policy22.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy22.proportionscreenedonadmission=1;
  policy22.proportionscreenedweekly=1;
  policy22.proportionscreenedondischarge=0;
  policy22.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken




  //  ===================================admission screening and weekly screening on a monday with an ideal test======================
//+ decolonisation of positives

  struct  interventionpolicy  policy23;

  policy23.screeninghighriskonly = false;
  policy23.admissionscreening.push_back(ideal1);
  policy23.weekdayscreening[0].push_back(ideal1);
  policy23.clinicalscreening.push_back(ideal1);
  //currently no weeklypostadmission or discharge screening
  policy23.preemptiveisolationforall=false;
  policy23.preemptiveisolationforhighrisk=false;
  policy23.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy23.isolateifpositive=false;
  policy23.blanketdecolforall=false;
  policy23.decolforhighrisk=false;
  policy23.decolifpositive=true;
  policy23.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy23.proportionscreenedonadmission=1;
  policy23.proportionscreenedweekly=1;
  policy23.proportionscreenedondischarge=0;
  policy23.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken

//  ===================================screening high risk only (with  an ideal test)======================
//===================================== + decolonisation for those who are positive ===========================================================

  struct  interventionpolicy  policy24;

  policy24.screeninghighriskonly = true;
  policy24.admissionscreening.push_back(ideal1);

  //NOT SURE THIS IS WORKING

  policy24.weekdayscreening[0].push_back(ideal1);
  policy24.clinicalscreening.push_back(ideal1);

  //currently no weeklypostadmission or discharge screening
  policy24.preemptiveisolationforall=false;
  policy24.preemptiveisolationforhighrisk=false;
  policy24.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy24.isolateifpositive=false;
  policy24.blanketdecolforall=false;
  policy24.decolforhighrisk=false;
  policy24.decolifpositive=true;
  policy24.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy24.proportionscreenedonadmission=1;
  policy24.proportionscreenedweekly=1;
  policy24.proportionscreenedondischarge=0;
  policy24.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken


//  ===================================no screening (just clinical swabs)======================
//===================================== just decolonisation for clinical cases ===========================================================

  struct  interventionpolicy  policy25;

  policy25.screeninghighriskonly = false;
  policy25.clinicalscreening.push_back(convculture1);

  policy25.preemptiveisolationforall=false;
  policy25.preemptiveisolationforhighrisk=false;
  policy25.numbernegativeuntilconsideredcleared=3; //i.e. 3 negative screens following a positive to confirm negative and
                                                        //change awareness state
  policy25.isolateifpositive=false;
  policy25.blanketdecolforall=false;
  policy25.decolforhighrisk=false;
  policy25.decolifpositive=true;
  policy25.timesteptoimplement=365;//i.e. at what timestep does the policy kick in
  policy25.proportionscreenedonadmission=1;
  policy25.proportionscreenedweekly=1;
  policy25.proportionscreenedondischarge=0;
  policy25.delaybeforeclinicalswab = 1;//when someone becomes infected a clinical swab will be taken - this is the time delay betweeen becoming I an dthis swab being taken




 switch(policynumber){
  case 1:
    policy=policy1;
    break;
  case 2:
    policy=policy2;
    break;
  case 3:
    policy=policy3;
    break;
  case 4:
    policy=policy4;
    break;
  case 5:
    policy=policy5;
    break;
  case 6:
    policy=policy6;
    break;
  case 7:
    policy=policy7;
    break;
  case 8:
    policy=policy8;
    break;
  case 9:
    policy=policy9;
    break;
  case 10:
    policy=policy10;
    break;
  case 11:
    policy=policy11;
    break;
  case 12:
    policy=policy12;
    break;
  case 13:
    policy=policy13;
    break;
  case 14:
    policy=policy14;
    break;
  case 15:
    policy=policy15;
    break;
  case 16:
    policy=policy16;
    break;
  case 17:
    policy=policy17;
    break;
  case 18:
    policy=policy18;
    break;
  case 19:
    policy=policy19;
    break;
  case 20:
    policy=policy20;
    break;
  case 21:
    policy=policy21;
    break;
  case 22:
    policy=policy22;
    break;
  case 23:
    policy=policy23;
    break;
  case 24:
    policy=policy24;
    break;
  case 25:
    policy=policy25;
    break;

  default:
    cout <<"Policy not specified so defauilting to policy 1";
     policy=policy1; //policy is the current active policy....this line sets default current policy to policy 1

  } //end switch

  policy.primary_isolation_effectiveness= primary_iso_effectiveness;
  policy.secondary_isolation_effectiveness= secondary_iso_effectiveness;
  policy.bodywash_pdc_effectiveness = bwash_pdc_effectiveness;
  policy.bodywash_pdi_effectiveness = bwash_pdi_effectiveness;
  policy.bodywash_progprob_effectiveness = bwash_progprob_effectiveness;
  policy.bodywash_IQ_effectiveness=0;
  policy.primary_isolation_capacity=get_isocap();

}

void patients::printpolicy(){
 //prints  the intervention policy
 cout<<"printpolicy not currently  implemented";


}




///////////////////////////method to read in files////////////////////////////////////////////

void read_in_files (char* Pdcfile, char* Pdifile, char* Progprobfile ) {//called at start of main i.e. just once for all simulations

cout<<"*****CALLINGFILES********"<< " \n";

/*ifstream re_ad_prob("daily readmission probs .txt");


   // read in header line from data file but don't use it
   char rubbish[256];
    re_ad_prob.getline(rubbish,256);


    //creates vector from read in file
    string y;
    while(getline(re_ad_prob, y)){

      re_Prob.push_back(atof(y.c_str()));

    }

    //cout<<"***Pdcfile: "<<Pdcfile<<" \n";
 // cout<<"***Pdifile: "<<Pdifile<<" \n";
  //cout<<"***Progprobfile: "<<Progprobfile<<" \n";

*/
////Readmission probabilities for SUS and COL



  ///////////////DISCHARGE AND DEATH PROBABILITIES FOR SUS AND INF
    ///////////// \ probabilities for susceptibles ///

    ifstream dis_prob_data_sus("daily discharge probs for the MRSA uninfected unadjusted.txt");


    // read in header line from data file but don't use it
    char rubbish[256];
    dis_prob_data_sus.getline(rubbish,256);


    //creates vector from read in file
    string x;
    while(getline(dis_prob_data_sus, x)){

      dis_prob_SUS.push_back(atof(x.c_str()));

    }



    /**
     //prints out the discaharge prob array
     for (unsigned int i = 0; i<dis_prob_SUS.size(); ++i){

     cout<<"\n Discharge probability array \n";
     cout<<"\n"<< dis_prob_SUS[i]<<"     printing our disharge probability array"<<"\n";

     }

    */

    /////////////daily death probabilities for susceptibles///////////////////////////////////////


      ifstream death_prob_data_sus("daily death probs for the MRSA uninfected unadjusted.txt");


      // read in header line from data file but don't use it            char rubbish[256];
      death_prob_data_sus.getline(rubbish,256);

      //creates vector from read in file
      string z;
      while(getline(death_prob_data_sus, z)){

 death_prob_SUS.push_back(atof(z.c_str()));


      }


      /**
       //prints out the death prob array
       for (unsigned int i = 0; i<death_prob_SUS.size(); ++i){

       cout<<"\n Death probability array \n";

       //cout<<death_prob_SUS[i]<<"printing our disharge probability array"<<"\n";

       }
      */


      ///////////// daily discharge probabilities for infecteds ///

        ifstream dis_prob_data_inf("daily discharge probs for the MRSA infected unadjusted.txt");
 //AT THE MINUTE THIS IS IDENTICAL TO TO THE DISCHARGE ARRAY FOR SUSCEPTIBLES


 // read in header line from data file but don't use it
 dis_prob_data_inf.getline(rubbish,256);

 //creates array from read in file
 string xx;
 while(getline(dis_prob_data_inf, xx)){

   dis_prob_INF.push_back(atof(xx.c_str()));

 }


 /**

 //prints out the discharge prob array
 for (unsigned int i = 0; i<dis_prob_INF.size(); ++i){

 cout<<"\n Discharge probability array for mrsa infected  \n";

 cout<<dis_prob_SUS[i]<<"... printing our disharge probability array"<<"\n";

 }
 */



        /////////////daily death probabilities for susceptibles///////////////////////////////////////

   ifstream death_prob_data_inf("daily death probs for the MRSA infected unadjusted.txt");


   // read in header line from data file but don't use it            char rubbish[256];
   death_prob_data_inf.getline(rubbish,256);

   //creates vector from read in file
   string yy;
   while(getline(death_prob_data_inf, yy)){

     death_prob_INF.push_back(atof(yy.c_str()));

   }


   /**
    //prints out the death prob array
    for (unsigned int i = 0; i<death_prob_INF.size(); ++i){

    cout<<"\n Death probability array for MRSA infected \n";

    cout<<death_prob_SUS[i]<<"...printing our death probability array"<<"\n";

    }
   */
   if(Pdcfile!=NULL){  //so there is a file to

     ifstream Pdc_data_file(Pdcfile);


     if (Pdc_data_file.bad()){
       cerr <<"Error: problem with "<<Pdcfile<<" \n";
       exit(8);
     }

     if(Pdc_data_file.is_open()){
       if(verbose)  cout <<"Pdc file is open" ;

       // read in header line from data file but don't use it
       Pdc_data_file.getline(rubbish,256);

       //creates vector from read in file
       string y;
       while(getline(Pdc_data_file, y)){
  //cout<<x<<"read in sens "<<endl;
  //cout<<sens_vector.size() <<"   size of sens_vector before push_back   "<<"\n";

  Pdc_vector.push_back(atof(y.c_str()));

  //cout<<Pdc_vector.size() <<"   size of Pdc_vector after push_back  "<<"\n";
       }

       //prints out the Pdc_vector
       if(verbose){
        for (unsigned int i = 0; i<Pdc_vector.size(); ++i){

  // if (DEFAULT_DEBUG==128){
   cout<<Pdc_vector[i]<<" Pdc  vector element"<<"\n";
  // }
        }
       }

     } else {// ! (Pdc_data_file.is_open()){
       //       cerr <<"Error: problem with Pdcfile "<<Pdcfile<<" File coult not be opened \n";
       cerr <<"Error: problem with Pdcfile "<<Pdcfile<<" File coult not be opened "<<" Pdc  value used is "<<Pdc<<" \n";
     }
   }

   if(Pdifile!=NULL){  //so there is a file to readin
     ifstream Pdi_data_file(Pdifile);

     if (Pdi_data_file.bad()){
       cerr <<"Error: problem with "<<Pdifile<<" \n";
       exit(8);
     }

     if( Pdi_data_file.is_open()){
       if(verbose) cout <<"Pdi file is open" ;


       // read in header line from data file but don't use it
       Pdi_data_file.getline(rubbish,256);

       //creates vector from read in file
       string y;
       while(getline(Pdi_data_file, y)){
  //cout<<x<<"read in sens "<<endl;
  //cout<<sens_vector.size() <<"   size of sens_vector before push_back   "<<"\n";

  Pdi_vector.push_back(atof(y.c_str()));

  //cout<<Pdi_vector.size() <<"   size of Pdi_vector after push_back  "<<"\n";
       }

       //prints out the Pdi_vector
              if(verbose){
           for (unsigned int i = 0; i<Pdi_vector.size(); ++i){

   //  if (DEFAULT_DEBUG==128){
   cout<<Pdi_vector[i]<<" Pdi from vector"<<"\n";
        }
       }
       // }
     } else {// ! (Pdi_data_file.is_open()){
       cerr <<"Error: problem with Pdifile "<<Pdifile<<" File coult not be opened "<<" Pdi  value used is "<<Pdi<<" \n";
     }

   }

   if(Progprobfile!=NULL){  //so there is a file to readin
     ifstream Progprob_data_file(Progprobfile);

     if (Progprob_data_file.bad()){
       cerr <<"Error: problem with "<<Progprobfile<<" \n";
       exit(8);
     }

     if(Progprob_data_file.is_open()){
       if(verbose) cout <<"Pdi file is open" ;

       // read in header line from data file but don't use it
       Progprob_data_file.getline(rubbish,256);

       //creates vector from read in file
       string y;
       while(getline(Progprob_data_file, y)){
  //cout<<x<<"read in sens "<<endl;
  //cout<<sens_vector.size() <<"   size of  before push_back   "<<"\n";

  Progprob_vector.push_back(atof(y.c_str()));

  //cout<<Progprob_vector.size() <<"   size of Progprob_vector after push_back  "<<"\n";
       }

       //prints out the Progprob_vector
       if(verbose){
  for (unsigned int i = 0; i<Progprob_vector.size(); ++i){

    //   if (DEFAULT_DEBUG==128){
    cout<<Progprob_vector[i]<<" Progprob from vector"<<"\n";
    //   }
  }
       }
     } else {// ! (Progprob_data_file.is_open()){
       cerr <<"Error: problem with Progprobfile "<<Progprobfile<<" File coult not be opened"<<" Progprob value used is "<<ProgProb<<" \n";
     }

   }

   /**


   ///////////////SENSITIVITY AND SPECIFICITY/////////////////////////////

   ////sensitivity
   ifstream sens_data_file("sensitivity_PCR.txt");

   if (sens_data_file.bad()){
   cerr <<"Error: problem with test_sens.dat \n";
   exit(8);
   }

   // read in header line from data file but don't use it
   sens_data_file.getline(rubbish,256);

   //creates vector from read in file
   string y;
   while(getline(sens_data_file, y)){
   //cout<<x<<"read in sens "<<endl;
   //cout<<sens_vector.size() <<"   size of sens_vector before push_back   "<<"\n";

   sens_vector.push_back(atof(y.c_str()));

   //cout<<sens_vector.size() <<"   size of sens_vector after push_back  "<<"\n";
   }

   //prints out the sens_vector
   for (unsigned int i = 0; i<sens_vector.size(); ++i){

   if (DEFAULT_DEBUG==128){
   cout<<sens_vector[i]<<" sens from vector"<<"\n";
   }
   }






   ////specificty
   ifstream spec_data_file("specificity_PCR.txt");

   if (spec_data_file.bad()){
   cerr <<"Error: problem with test_spec.dat \n";
   exit(8);
   }

   // read in header line from data file but don't use it
   spec_data_file.getline(rubbish,256);

   //creates vector from read in file

   string u;
   while(getline(spec_data_file, u)){
   //cout<<u<<"read in spec"<<endl;
   //cout<<"in while loop"<<"\n";
   //cout<<spec_vector.size() <<"   size of spec_vector before push_back   "<<"\n";

   spec_vector.push_back(atof(u.c_str()));

   //cout<<spec_vector.size() <<"   size of spec_vector after push_back  "<<"\n";
   }

   //prints out the SandS_vector
   for (unsigned int i = 0; i<spec_vector.size(); ++i){
   if (DEFAULT_DEBUG==128){
   cout<<spec_vector[i]<<"  spec from vector"<<"\n";
   }
   }

   */

   /**

   ///////////////DEMOGRAPHICS/////////////////////////////
   //not using this at the moment
   ifstream age_data_file("test_age.txt");

   if (age_data_file.bad()){
   cerr <<"Error: problem with test_age.dat \n";
   exit(8);
   }

   // read in header line from data file but don't use it
   age_data_file.getline(rubbish,256);


   ////creates vector from read in file/
   string w;
   while(getline(age_data_file, w)){
   //cout<<"in while loop"<<"\n";
   //cout<<age_vector.size() <<"   size of age_vector before push_back   "<<"\n";
   //cout<<atoi(w.c_str());
   age_vector.push_back(atoi(w.c_str()));
   //cout<<age_vector.size() <<"   size of age_vector after push_back  "<<"\n";
   }

   //prints out the age_vector
   for (unsigned int i = 0; i<age_vector.size(); ++i){
   if (DEFAULT_DEBUG==128){
   cout<<age_vector[i]<<"age from vector"<<"\n";
   }
   }





   //counting numbers in age groups
   for (unsigned long int q=0; q<age_vector.size(); q++){

   //cout<<int(age_vector[q])<< "age vector [q] "<<"\n";
   if (int(age_vector[q])<=19){
   a++;}
   //cout<<a<<" number in group a"<<"\n";
   if (int(age_vector[q])>=20 && int(age_vector[q])<=39){
   b++;}
   //cout<<b<<" number in group b"<<"\n";
   if (int(age_vector[q])>=40 && int(age_vector[q])<=59){
   c++;}
   //cout<<c<<" number in group c"<<"\n";
   if (int(age_vector[q])>=60 && int(age_vector[q])<=79){
   d++;}
   //cout<<d<<" number in group d"<<"\n";
   if (80<=int(age_vector[q])){
   e++;}
   //cout<<e<<" number in group e"<<"\n";

   }

   double totalnopeople=a+b+c+d+e;
   //cout<<totalnopeople<<"   total no people   "<<a<<"  a"<<"\n";

   proba = a/totalnopeople;
   //cout<<proba<<"proba"<<"\n";
   probb = proba+b/totalnopeople;
   //cout<<probb<<"probb"<<"\n";
   probc = probb+c/totalnopeople;
   //cout<<probc<<"probc"<<"\n";
   probd = probc+d/totalnopeople;
   //cout<<probd<<"probd"<<"\n";
   probe = 1;  //e/totalnopeople;
   //cout<<probe<<"probe"<<"\n";


   //LoS_data_file.close();
   sens_data_file.close();
   spec_data_file.close();
   age_data_file.close();

   */

   dis_prob_data_sus.close();
   death_prob_data_sus.close();
   dis_prob_data_inf.close();
   death_prob_data_inf.close();
   //Pdc_data_file.close();
   //Pdi_data_file.close();
   //Progprob_data_file.close();

}//end of read_in_files




//////////////////////////method to set parameters to contents of file (for each simulation run)/////////////////////

    void read_in_parameters(gsl_rng * rng){//called within for loop in main - therefore reset for each simulation run
        //cout<<"in read in parameters"<<"\n";

            float urn;
            urn=gsl_rng_uniform(rng);

            if (DEFAULT_DEBUG==128){
            cout<<sens_vector.size()<<"sens_vector.size()"<<"\n";
            }


            int chosen_SandS = int(urn * (sens_vector.size()+1));


            if (DEFAULT_DEBUG==128){
            cout<<chosen_SandS<<"chosen_SandS"<<"\n";
            }

//            sensitivity = sens_vector[chosen_SandS];
//           specificity = spec_vector[chosen_SandS];

            //cout<<sensitivity<<"  sensitivity"<<"\n";
            //cout<<specificity<<"  specificity"<<"\n";

    }


     bool patients::isthereweeklypostadmissionscreening(void){ //true if there is, false if there isn't
       if(policy.weeklypostadmission.begin()==policy.weeklypostadmission.end()){
         return(false); //since no weekly post admsisin screening
       } else {
         return(true);
       }
    };

    int patients::get_weekday_number(void){
    //return day of week: 0 for monday, 1 for tuesday etc (accuonting for timesteps per day).

        int  daynum; // accounting for timesteps perday
        int  dayofweeknum;

        daynum=time/timestepsperday;
        //cout<<"daynum="<<daynum;
        dayofweeknum=daynum % 7;
        return(dayofweeknum);
    };

     bool patients::is_it_first_time_step_of_day(void){
     // returns true if yes false if not
       int daynum1;
       int daynum2 ;
       daynum2=time/timestepsperday;
       daynum1=(time-1)/timestepsperday;
       if(daynum1!=daynum2) {
            //cout <<"first time step of day\n";
            return(true);
        } else {
            return(false);
        }
     };

     //set up output files



 void patients::run_n_simulations(const int stoptime, const int numinitialcases, gsl_rng *rng, const long int numsims){





    //  stepresults = new ICUstate*[numsims];


  /*for(int simnum=0; simnum<numsims;simnum++)
	{
		// allocate realization.
		stepresults[simnum] = new ICUstate[stoptime];
		memset(stepresults[simnum], 0, sizeof(ICUstate)*stoptime);

		//initialise population in each realization.
		stepresults[simnum][0].SinICU = 0;
		stepresults[simnum][0].CinICU = 0;
		stepresults[simnum][0].IinICU = 0;
		stepresults[simnum][0].ISOinICU = 0;
		stepresults[simnum][0].time = 0;
	}
*/

//  char results_sample[]= "results_sample.txt";
//  std::ofstream resultsStream(results_sample);

 char results_timesteps[]= "results_timesteps.txt";
 std::ofstream resultsStream(results_timesteps);

 char prev_file[]="prevalence_over_time.txt";
 std::ofstream prevStream(prev_file);

 char homeduration_file[]="homeduration.txt";
 std::ofstream durStream(homeduration_file);

//use these to send detailed output directly to txt files (if required)
//ofstream prev_file("prevalence_over_time.txt",  ios::app );
//ofstream bed_days_file("bed_days_in_ICU.txt",  ios::app );
//cout<<"let's go.........................................................\n\n\n\n";

  double meanS =0;
  double meanC =0;
  double meanI =0;
  double meanISO =0;

  SinICU=0;
  CinICU=0;
  IinICU=0;
  ISOinICU=0;
  HIGHRISKinICU=0;
  noICU=0;


  float meansusbedday =0;
  float meancolbedday =0;
  float meaninfbedday=0;
  float meanisobedday=0;
  float meandecolcount=0;
  float meanno_ad_screens=0;
  float meanno_wkly_screens=0;
  float meanno_clin_screens=0;
  float meanno_pos_screens=0;
  float meanno_neg_screens=0;
  float meanno_pos_CC_screens=0;
  float meanno_neg_CC_screens=0;
  float meanno_pos_CA_screens=0;
  float meanno_neg_CA_screens=0;
  float meanno_pos_CA_early_screens=0;
  float meanno_neg_CA_early_screens=0;
  float meanno_PCR_screens=0;
  float meancumulative_StoC=0;
  float meancumulative_StoI =0;
  float meancumulative_CtoI =0;
  float meancumulative_admissions =0;
  float meancumulative_discharges=0;
  float meancumulative_deaths=0;
        float meancumulative_appisodays =0;
  float meancumulative_inappisodays=0;
  float meancumulative_unisodays=0;
  float meancumulative_good_bd=0;
  float meanbeddaycosts=0;
  float meanisolationcosts=0;
  float meandecolcosts =0;
  float meanswabbingcosts=0;
  float meanscreeningcosts=0;
  float meantreatmentcosts=0;
  float meantotalcosts=0;
  float meanhealthbenefitsinICU=0;


        //start off with  colonized individuals in the ICU
     //number of these is given by numinitial cases (which can be set to 0 in main
     //in which case instead of having colonized intially in the ward they are imported into it in the admission routine

   //  unsigned long int persontoinfect; //identifier for the intial people to infect
//  struct events startevent = {COLONIZE, NULL}; //infect event character codes. Member order: event type, pointer to patient

       // don't need this since there is some probabity people will be colonized on admission, so no need to start with any cases
//  for(int i=0; i<numinitialcases; ++i){
//
//   persontoinfect=int((DEFAULT_ICU_SIZE-0.0000001)*gsl_rng_uniform(rng));//randomly choose a person from the popn to infect
//   startevent.patient_ptr=ICU_patients[persontoinfect].ICU_patient_ptr;//make the pointer=the address of the randomly chosen person
//
//   // infect these selected patients at time 0
//   time=0;
//   future_infection_events[time].push_back(startevent);//put the event onto the end of the list of infection events to do
//  }

  if (DEFAULT_DEBUG==103){
    resultsStream<<"Sim "<<" \t"<<"Day "<<" \t"<<" ICUsize "<<" \t"<<" Sus "<<" \t"<<" Col "<<" \t"<<" Inf "<<" \t"<<" Iso "<<" \t"<<" High Risk "<<" \t"<<" Sbd "<<" \t"<<" Cbd "<<" \t"<<" Ibd "<<" \t"<<" Isobd "<<" \t"<<" Decolcount "<<" \t"<<" no_ad_screens "<<" \t"<<" no weekly screens "<<" \t"<<" no clinical screens "<<" \t"<<" no positive screens "<<" \t"<<" no neg screens "<<" \n";

  }

  //prevStream<<"numsim"<<","<<"time"<<","<<"homeduration"<<","<<"numbercolonised"<<","<<"numberinfected"<<",\n";

  for(int simnum=0; simnum<numsims; ++simnum){



    //  ICUstate* currentrun = stepresults[simnum];

     // reset_all_counts();

   //need to empty all queues here as there will be some future events left in the queue from previous simulations which will slow things down
   //increasingly and waste memory as numsims increases
   clear_future_events();
   reset_all_counts();

   reset_all_counts();
   // cout<<" stoptime    is "<< stoptime<< "\n";
   //a simulation

        unsigned long int pick_patient;
        float rannum;
        rannum = gsl_rng_uniform(rng);//uniform random number: includes 0 but excludes 1
        pick_patient = int (DEFAULT_ICU_SIZE* rannum);



    while(hos_pop[pick_patient].ICU_no==99){ //i.e. keep trying to pick a patient until we find one not in ICU (so ICU_no==99)
        rannum = gsl_rng_uniform(rng);//uniform random number: includes 0 but excludes 1
        pick_patient = int (popsize * rannum);

 }


  int sample =pick_patient;

   for(time=0; time<stoptime; ++time){ //time specifies the time step (not necessarily the same as actual time is there me be one or more timesteps per day)


         if (DEFAULT_DEBUG==8){

         if (time==policy.timesteptoimplement){

             cout<<"*****************SCREENING AND INTERVENTION BEGINS************************"<<" \n";

         }
         }

   if (DEFAULT_DEBUG==100){




      resultsStream<<  "\n"<< hos_pop[sample].patientid<<" SAMPLE-PATIENT ID"<<"\n";
      resultsStream <<    hos_pop[sample].no_hos<< " SAMPLE-Number Hospitilisations"<<"\n";

        if (hos_pop[sample].disease_state==COLONIZED) {

       resultsStream<< " COLONIZED"<<"\n";
        }

        if (hos_pop[sample].disease_state==SUSCEPTIBLE){

       resultsStream<< " SUSCEPTIBLE"<<"\n";
        }

        if (hos_pop[sample].disease_state==INFECTED){

       resultsStream<< " INFECTED"<<"\n";
        }

        if (hos_pop[sample].specialty==GM){

            resultsStream<< " SAMPLE-Location GM"<<"\n";

        }

        if (hos_pop[sample].specialty==COM){

           resultsStream<< " SAMPLE-Location Community"<<"\n";

       }

        if (hos_pop[sample].specialty==ICU){

            resultsStream<< " SAMPLE-Location ICU"<<"\n";

        }



         if (hos_pop[sample].risk_group==HIGH){
        resultsStream<< " SAMPLE-HIGH Risk group"<<"\n"<<"\n";

         }


         if (hos_pop[sample].risk_group==LOW){
         resultsStream<< " SAMPLE-LOW Risk group"<<"\n"<<"\n";

         }
      }



susbedday=0;
colbedday=0;
infbedday=0;
isobedday=0;



   /* int GENCOL=0;

    for (int j =0; j< DEFAULT_POPSIZE; ++j)

    {
        if (hos_pop[j].disease_state==COLONIZED){

        ++GENCOL;}
    }
    */

    SinICU=0;
    CinICU=0;
    IinICU=0;
    ISOinICU=0;
    HIGHRISKinICU=0;
    noICU=0;

    //cout<<"Time is "<<time<<" Day is "<<get_weekday_number()<<"\n";
                bool firsttsofday=is_it_first_time_step_of_day(); //true if this is the first time step of the day

                long int burn_in_time = 0; //during the burn-in period (specified here) results aren't stored - allow system to reach equilibrium

    if (time==burn_in_time){ //once the burn-in period is up reset all output variables and start counting again from here.
   //  reset_all_counts();
                }

    process_infections(rng); //see if new patients should become infected/colonized, and add to schedule
    //print_events();

    if(firsttsofday)  process_movements(rng); //i.e. only see if patients discharge or die on first timestep of day

    if (DEFAULT_DEBUG==101)
    {       print_ICU_data();
    }  // This prints ICU state on the current day
    // update current and cumulative appropriate and inappropriate isolation days and unisolated days
    appisodays=0;
    inappisodays=0;
    unisodays=0;
    app_or_inapp_isolation_days(); //calculates current appropriate and inappripriate isolation days and updates cumulative totoals

    bed_day_counts();

    LOS_counts();

    //perform all day of week screens (i.e. regular monday screens etc)
    //note that this code is quite inefficient...so if speed becomes an issue can optimise this


  if (time>policy.timesteptoimplement){

    for (int j =0; j< DEFAULT_ICU_SIZE; ++j) {
     if (ICU_patients[j].ICU_patient_ptr!=NULL){ //i.e. if there is someone in the bed
      scheduled_dayofweek_screening(ICU_patients[j].ICU_patient_ptr, rng);
     }
    }

  }

    perform_events(rng);





    for (int j =0; j< DEFAULT_ICU_SIZE; ++j){//20 should be size of ICU_patients

                    noICU++;



//                if (ICU_patients[j].ICU_patient_ptr->disease_state==SUSCEPTIBLE|ICU_patients[j].ICU_patient_ptr->disease_state==RECOVEREDFROMCOL|ICU_patients[j].ICU_patient_ptr->disease_state==RECOVEREDFROMINF){
                if (ICU_patients[j].ICU_patient_ptr->disease_state==SUSCEPTIBLE){

                    ++SinICU;

     // currentrun[time].SinICU=SinICU;

                    }

                    if (ICU_patients[j].ICU_patient_ptr->disease_state==COLONIZED){

                        ++CinICU;

     // currentrun[time].CinICU=CinICU;
                    }

     //                if (ICU_patients[j].ICU_patient_ptr->disease_state==INFECTED|ICU_patients[j].ICU_patient_ptr->disease_state==INFECTEDFROMCOL){
                    if (ICU_patients[j].ICU_patient_ptr->disease_state==INFECTED){

                        ++IinICU;
    //  currentrun[time].IinICU=IinICU;
                    }

                    if (ICU_patients[j].ICU_patient_ptr->isolation==true){
                        ++ISOinICU;
      //currentrun[time].ISOinICU=ISOinICU++;
                    }


             if (ICU_patients[j].ICU_patient_ptr->no_hos >1){


                  ++ HIGHRISKinICU;

             }

//cout<<"  "<<time<<"  "<<SinICU<<"  "<<CinICU<<"  "<<IinICU<<"  "<<ISOinICU<<"\n ";




    /**
     CUMULATIVE OUTPUTS

     */



     }





if (DEFAULT_DEBUG==103){
    //resultsStream<<"Day "<<time<<" Sus "<<S<<" Col "<<C<<" Inf "<<I<<" iso "<<ISO<<" Sbd "<<susbedday<<" Cbd "<<colbedday<<" Ibd "<<infbedday<<" Isobd "<<isobedday<<" Decolcount "<<decolcount<<" no_ad_screens "<<no_ad_screens<<" no weekly screens "<<no_wkly_screens<<" no clinical screens "<<no_clin_screens<<" no positive screens "<<no_pos_screens<<" no neg screens "<<no_neg_screens<<" \n";

    resultsStream<<simnum<<" \t"<<time<<" \t"<<noICU<<" \t"<<SinICU<<" \t"<<CinICU<<" \t"<<IinICU<<" \t"<<ISO<<" \t"<<HIGHRISKinICU<<" \t"<<susbedday<<" \t"<<colbedday<<" \t"<<infbedday<<" \t"<<isobedday<<" \t"<<decolcount<<" \t"<<no_ad_screens<<" \t"<<no_wkly_screens<<" \t"<<no_clin_screens<<" \t"<<no_pos_screens<<" \t"<<no_neg_screens<<" \n";



    }

    if (DEFAULT_DEBUG==104){
    //resultsStream<<"Day "<<time<<" Sus "<<S<<" Col "<<C<<" Inf "<<I<<" iso "<<ISO<<" Sbd "<<susbedday<<" Cbd "<<colbedday<<" Ibd "<<infbedday<<" Isobd "<<isobedday<<" Decolcount "<<decolcount<<" no_ad_screens "<<no_ad_screens<<" no weekly screens "<<no_wkly_screens<<" no clinical screens "<<no_clin_screens<<" no positive screens "<<no_pos_screens<<" no neg screens "<<no_neg_screens<<" \n";

 //   resultsStream<<simnum<<" \t"<<time<<" \t"<<CtoI<<" \t"<<Conadmissionif (DEFAULT_DEBUG==103){
    //resultsStream<<"Day "<<time<<" Sus "<<S<<" Col "<<C<<" Inf "<<I<<" iso "<<ISO<<" Sbd "<<susbedday<<" Cbd "<<colbedday<<" Ibd "<<infbedday<<" Isobd "<<isobedday<<" Decolcount "<<decolcount<<" no_ad_screens "<<no_ad_screens<<" no weekly screens "<<no_wkly_screens<<" no clinical screens "<<no_clin_screens<<" no positive screens "<<no_pos_screens<<" no neg screens "<<no_neg_screens<<" \n";

 //   resultsStream<<simnum<<" \t"<<time<<" \t"<<noICU<<" \t"<<SinICU<<" \t"<<CinICU<<" \t"<<IinICU<<" \t"<<ISO<<" \t"<<HIGHRISKinICU<<" \t"<<susbedday<<" \t"<<colbedday<<" \t"<<infbedday<<" \t"<<isobedday<<" \t"<<decolcount<<" \t"<<no_ad_screens<<" \t"<<no_wkly_screens<<" \t"<<no_clin_screens<<" \t"<<no_pos_screens<<" \t"<<no_neg_screens<<" \n";if (DEFAULT_DEBUG==103){
    //resultsStream<<"Day "<<time<<" Sus "<<S<<" Col "<<C<<" Inf "<<I<<" iso "<<ISO<<" Sbd "<<susbedday<<" Cbd "<<colbedday<<" Ibd "<<infbedday<<" Isobd "<<isobedday<<" Decolcount "<<decolcount<<" no_ad_screens "<<no_ad_screens<<" no weekly screens "<<no_wkly_screens<<" no clinical screens "<<no_clin_screens<<" no positive screens "<<no_pos_screens<<" no neg screens "<<no_neg_screens<<" \n";

 //   resultsStream<<simnum<<" \t"<<time<<" \t"<<noICU<<" \t"<<SinICU<<" \t"<<CinICU<<" \t"<<IinICU<<" \t"<<ISO<<" \t"<<HIGHRISKinICU<<" \t"<<susbedday<<" \t"<<colbedday<<" \t"<<infbedday<<" \t"<<isobedday<<" \t"<<decolcount<<" \t"<<no_ad_screens<<" \t"<<no_wkly_screens<<" \t"<<no_clin_screens<<" \t"<<no_pos_screens<<" \t"<<no_neg_screens<<" \n";



    }

    }// end for time




    //COSTS////////////////////////////////////////////


     if (DEFAULT_DEBUG==8){
            cout<<ISO<<"**VALUE of ISO count at end of simulation loop"<<"\n";
        }





   beddaycosts = getbeddaycosts();
   isolationcosts = getisolationcosts();
   decolcosts=getdecolcosts();
   swabbingcosts = getswabbingcosts();
   screeningcosts = getscreeningcosts();
   treatmentcosts = gettreatmentcosts();

   totalcosts = beddaycosts+ isolationcosts + decolcosts + swabbingcosts + screeningcosts+treatmentcosts;
   healthbenefitsinICU=gethealthbenefitsinICU();
            /////////////////////////////////////HEALTH BENEFITS////////////////////////////////////////////


   /**  -----------------  OUTPUT -------------------*/

   if(BATCHMODE>0){
//julie debug   if(DEFAULT_DEBUG>=1){ //then output summary details for every sim
     std::cout<<susbedday<<" "<<colbedday<<" "<<infbedday<<" "<<isobedday<<" "<<decolcount<<" "<<no_ad_screens<<" "<<no_wkly_screens<<" "<<no_clin_screens<<" "<<no_pos_screens<<" "<<no_neg_screens;
     std::cout<<" "<<no_pos_CC_screens<<" "<<no_neg_CC_screens<<" "<<no_pos_CA_screens<<" "<<no_neg_CA_screens;
     std::cout<<" "<<no_pos_CA_early_screens<<" "<<no_neg_CA_early_screens<<" "<<no_PCR_screens;
     std::cout<<" "<<cumulative_StoC<<" "<<cumulative_StoI<<" "<<cumulative_CtoI <<" "<<cumulative_admissions<<" "<<cumulative_discharges<<" "<< cumulative_deaths;
     std::cout<<"  "<<cumulative_appisodays <<" "<<cumulative_inappisodays<<"  "<<cumulative_unisodays <<"  "<<cumulative_good_bd;
    std::cout<<"  "<<beddaycosts<<"  "<<isolationcosts<<"  "<<decolcosts<<"  "<<swabbingcosts<<"  "<<screeningcosts<<"  "<<treatmentcosts<<"  "<<totalcoldischarged <<" "<<totalcosts;
    std::cout<<"  "<<healthbenefitsinICU<<" \n";//
 //   }
   }


    meansusbedday+=susbedday;
    meancolbedday+=colbedday;
    meaninfbedday+=infbedday;
    meanisobedday+=isobedday;
    meandecolcount+=decolcount;
    meanno_ad_screens+=no_ad_screens;
    meanno_wkly_screens+=no_wkly_screens;
    meanno_clin_screens+=no_clin_screens;
    meanno_pos_screens+=no_pos_screens;
    meanno_neg_screens+=no_neg_screens;
    meanno_pos_CC_screens+=no_pos_CC_screens;
    meanno_neg_CC_screens+=no_neg_CC_screens;
    meanno_pos_CA_screens+=no_pos_CA_screens;
    meanno_neg_CA_screens+=no_neg_CA_screens;
    meanno_pos_CA_early_screens+=no_pos_CA_early_screens;
    meanno_neg_CA_early_screens+=no_neg_CA_early_screens;
    meanno_PCR_screens+=no_PCR_screens;
    meancumulative_StoC+=cumulative_StoC;
    meancumulative_StoI+=cumulative_StoI;
    meancumulative_CtoI+=cumulative_CtoI;
    meancumulative_admissions+=cumulative_admissions;
    meancumulative_discharges+=cumulative_discharges;
    meancumulative_deaths+=cumulative_deaths;
    meancumulative_appisodays+=cumulative_appisodays;
    meancumulative_inappisodays+=cumulative_inappisodays;
    meancumulative_unisodays+=cumulative_unisodays;
    meancumulative_good_bd+=cumulative_good_bd;
    meanbeddaycosts+=beddaycosts;
    meanisolationcosts+=isolationcosts;
    meandecolcosts +=decolcosts;
    meanswabbingcosts+=swabbingcosts;
    meanscreeningcosts+=screeningcosts;
    meantreatmentcosts+=treatmentcosts;
    meantotalcosts+=totalcosts;
    meanhealthbenefitsinICU+=healthbenefitsinICU;



            //std::cout<<"\nFinal data: \n";
       //
       //     bed_days_file.close();
//  prev_file.close();



} //end for simnum

if (DEFAULT_DEBUG==108){

 for(time=0; time<stoptime; ++time){

        meanS/=numsims;
        meanC/=numsims;
        meanI/=numsims;
        meanISO/=numsims;

        resultsStream<<time<<" \t"<<meanS<<" \t"<<meanC<<" \t"<<meanI<<" \t"<<meanISO<<" \n";

 }
}


//for(time=0; time<stoptime; ++time)
//	{
		//double prevC = 0;

//	for(int simnum=0; simnum<numsims; ++simnum)

	//	{
	//		prevC=stepresults[simnum][time].CinICU;


	//	}


		//prevC /= numsims;
		//prev_file << stepresults[3][time].time << "\t" << stepresults[3][time].CinICU  << "\n";


//	}


        meansusbedday/=numsims;
        meancolbedday/=numsims;
        meaninfbedday/=numsims;
        meanisobedday/=numsims;
        meandecolcount/=numsims;
        meanno_ad_screens/=numsims;
        meanno_wkly_screens/=numsims;
        meanno_clin_screens/=numsims;
        meanno_pos_screens/=numsims;
        meanno_neg_screens/=numsims;
        meanno_pos_CC_screens/=numsims;
  meanno_neg_CC_screens/=numsims;
  meanno_pos_CA_screens/=numsims;
  meanno_neg_CA_screens/=numsims;
  meanno_pos_CA_early_screens/=numsims;
  meanno_neg_CA_early_screens/=numsims;
  meanno_PCR_screens/=numsims;
  meancumulative_StoC/=numsims;
  meancumulative_StoI/=numsims;
  meancumulative_CtoI/=numsims;
  meancumulative_admissions/=numsims;
  meancumulative_discharges/=numsims;
  meancumulative_deaths/=numsims;
  meancumulative_appisodays/=numsims;
  meancumulative_inappisodays/=numsims;
  meancumulative_unisodays/=numsims;
  meancumulative_good_bd/=numsims;
  meanbeddaycosts/=numsims;
  meanisolationcosts/=numsims;
  meandecolcosts/=numsims;
  meanswabbingcosts/=numsims;
  meanscreeningcosts/=numsims;
  meantreatmentcosts/=numsims;

  meantotalcosts/=numsims;
//;
  meanhealthbenefitsinICU/=numsims;



 //resultsStream<<meansusbedday<<" "<<meancolbedday<<" "<<meaninfbedday<<" "<<meanisobedday<<" "<<meandecolcount<<" "<<meanno_ad_screens<<" "<<meanno_wkly_screens<<" "<<meanno_clin_screens<<" "<<meanno_pos_screens<<" "<<meanno_neg_screens;
 std::cout<<" "<<meanno_pos_CC_screens<<" "<<meanno_neg_CC_screens<<" "<<meanno_pos_CA_screens<<" "<<meanno_neg_CA_screens;
 std::cout<<" "<<meanno_pos_CA_early_screens<<" "<<meanno_neg_CA_early_screens<<" "<<meanno_PCR_screens;
 std::cout<<" "<<meancumulative_StoC<<" "<<meancumulative_StoI<<" "<<meancumulative_CtoI <<" "<<meancumulative_admissions<<" "<<meancumulative_discharges<<" "<< meancumulative_deaths;
 std::cout<<" "<<meancumulative_appisodays <<" "<<meancumulative_inappisodays<<"  "<<meancumulative_unisodays <<"  "<<meancumulative_good_bd;

//cout<<"pdc is "<<Pdc<<"\n";
 std::cout<<"  "<<meanbeddaycosts<<"  "<<meanisolationcosts<<"  "<<meandecolcosts<<"  "<<meanswabbingcosts<<"  "<<meanscreeningcosts<<"  "<<meantreatmentcosts<<"  "<<meantotalcosts;
 std::cout<<"  "<<meanhealthbenefitsinICU<<" \n";//



};//end of run_n_simulations

void patients::reset_all_counts(){

///run though and reset hos_pop-this seems to be neccessary to make each simulation begin afresh
for(unsigned long int i=0; i<popsize; ++i){


        hos_pop[i].disease_state=  SUSCEPTIBLE;
        hos_pop[i].awareness_state=UNKNOWN;
        hos_pop[i].no_hos=0; //count if ever in hospital(ICU)
        hos_pop[i].risk_group=LOW;

    }

S=popsize;//-100;
 C=0;//90
 I=0;//10
 CI=0;
 ISO=0;
 SECISO=0;
 D=0;


SinICU=0;
CinICU=0;
IinICU=0;
ISOinICU=0;
noICU=0;
HIGHRISKinICU=0;


 susbedday=0; colbedday=0; infbedday=0; isobedday=0; decolcount=0; no_ad_screens=0; no_wkly_screens=0; no_clin_screens=0; no_pos_screens=0; no_neg_screens=0;
 no_pos_CC_screens = 0; no_neg_CC_screens=0; no_pos_CA_screens=0; no_neg_CA_screens=0; no_pos_CA_early_screens=0; no_neg_CA_early_screens=0;
 no_PCR_screens = 0;
 cumulative_StoC=0; cumulative_StoI=0;
 cumulative_CtoI =0; cumulative_admissions=0; cumulative_discharges=0; cumulative_deaths=0;
 cumulative_appisodays=0; cumulative_inappisodays=0;cumulative_unisodays=0 ; cumulative_good_bd = 0; beddaycosts=0; isolationcosts=0; decolcosts=0; swabbingcosts=0; screeningcosts=0; treatmentcosts=0; totalcosts=0 ;
    healthbenefitsinICU=0;

}

/**
NOT USING THIS ANYMORE BECAUSE CALCULATING DISCHARGE AND DEATH PROBS DAILY RATHER THAN ASSIGNING A PRESET LOS ON ADMISSION
////////////////////////////method to randomly chose a LoS from the distribution (vector)/////////////////////////////

    int patients::pick_LoS(gsl_rng * rng){
        float urn;
        urn=gsl_rng_uniform(rng);

            //cout<<LoS_vector.size()<<"LoS_vector.size()"<<"\n";

            int chosen_LoS = (urn * LoS_vector.size());
            LoS = LoS_vector[chosen_LoS];
            //cout<<LoS<<"LoS"<<"\n";
            return LoS;

    }


*/





/**

    //int patients::pick_add_LoS(gsl_rng * rng){
//picking additional length of stay from a normal distribution (bayersman paper)
    int patients::pick_add_LoS(gsl_rng * rng){
*/



/**
        int pre_transformed_add_LoS;
        int sigma = 5;
        int mean = 5;
        pre_transformed_add_LoS = gsl_ran_gaussian (rng, sigma);

        int add_LoS;
        add_LoS = mean+pre_transformed_add_LoS;

//ensures that cannot take time off stay - thsi needs to be dealt with so that this is possible.
            if (add_LoS<=0){
                add_LoS=pick_add_LoS(rng);
            }
            if (add_LoS>0){

              cout<<add_LoS<<"add_LoS returned"<<"\n";

              return add_LoS;

            }
*/





/**
        int add_LoS;
        add_LoS=5;
        return add_LoS;
*/




/**
        //cout<<"pick add Los called"<<"\n";
        float urn;
        urn=gsl_rng_uniform(rng);



            int chosen_add_LoS = (urn * add_LoS_vector.size());
            add_LoS = add_LoS_vector[chosen_add_LoS];

            return add_LoS;
*/


///////////////////////////method to randomly chose a LoS from the distribution (vector)/////////////////////////////


//THIS ISNT CALLED AT THE MOMENT - NOT DEALING WITH AGE
    int patients::calculate_characteristics(gsl_rng * rng){

        //demographics of admission population



            float ran;
            ran=gsl_rng_uniform(rng);

            //cout<<ran<<"the random number"<<"\n";
            //cout<<proba<< "prob a"<<"\n";

            if (ran<proba){
                age = 19;//means they are in age group 0 -19
                }
                else {
                    if (ran<probb){
                        age = 2039;//age group 20-39
                    }
                        else {
                            if (ran<probc){
                                age = 4059;//age group 40-59
                            }
                                else {
                                    if (ran<probd){
                                    age = 6079;//age group 60-79
                                    }
                                    else{
                                        if (ran<probe){
                                        age = 80;//age group 80+
                                        }
                                    }
                                }
                        }
                }

            //cout<<age<< "age "<<"\n";
            return age;


    }



///////////////////////////////////////IQ calculations/////////////////////////////////


////////////baseline IQ
//this selects IQ from dist and returns it
//all IQs are at the minute set to 1 though
    double get_IQ(gsl_rng * rng){

        double baseline_IQ;
        double sigma = 1;
   //     double zeta = 0.5;
        baseline_IQ = gsl_ran_gaussian (rng, sigma);


        return baseline_IQ;


    }


    ////effect on IQ due to ISOLATION


// /**
//     double effect_of_ISO(gsl_rng * rng){
     //effect of isolation in reducing transmission from a colonized/infected isolated patient
     //Isolation may also reduce the risk of acquisition of MRSA amongst susceptible isolated patients
     // however this is not currently allowed for.

//         double ISO_effect;
//         ISO_effect=0.35;
//           ISO_effect= ISOEFFECT;
//         return ISO_effect;
//     }


// */

        double get_effect_of_ISO(gsl_rng * rng){

            double pre_transformed_ISO_effect;
            double sigma = effect_of_ISO_SD;
            double mean = effect_of_ISO;
            pre_transformed_ISO_effect = gsl_ran_gaussian (rng, sigma);

            double ISO_effect;
            ISO_effect = mean+pre_transformed_ISO_effect;


            //cout<<ISO_effect<<" ISO_effect "<<"\n";

            return ISO_effect;


        }


        double get_pdc_effect_of_bodywash(gsl_rng * rng){

            double pre_transformed_bodywash_pdc_effect;
            double sigma = DecolEffectonPdc_SD;
            double mean = DecolEffectonPdc_MEAN;
            pre_transformed_bodywash_pdc_effect = gsl_ran_gaussian (rng, sigma);

            double bodywash_pdc_effect;
            bodywash_pdc_effect = mean+pre_transformed_bodywash_pdc_effect;



            return bodywash_pdc_effect;


        }

        double get_pdi_effect_of_bodywash(gsl_rng * rng){

            double pre_transformed_bodywash_pdi_effect;
            double sigma = DecolEffectonPdi_SD;
            double mean = DecolEffectonPdi_MEAN;
            pre_transformed_bodywash_pdi_effect = gsl_ran_gaussian (rng, sigma);

            double bodywash_pdi_effect;
            bodywash_pdi_effect = mean+pre_transformed_bodywash_pdi_effect;



            return bodywash_pdi_effect;


        }



         double get_IQ_effect_of_bodywash(gsl_rng * rng){

            double pre_transformed_bodywash_IQ_effect;
            double sigma = DecolEffectonIQ_SD;
            double mean = DecolEffectonIQ_MEAN;
            pre_transformed_bodywash_IQ_effect = gsl_ran_gaussian (rng, sigma);

            double bodywash_IQ_effect;
            bodywash_IQ_effect = mean+pre_transformed_bodywash_IQ_effect;



            return bodywash_IQ_effect;


        }


        double get_progprob_effect_of_bodywash(gsl_rng * rng){

            double pre_transformed_bodywash_progprob_effect;
            double sigma = DecolEffectonProgProb_SD;
            double mean = DecolEffectonProgProb_MEAN;
            pre_transformed_bodywash_progprob_effect = gsl_ran_gaussian (rng, sigma);

            double bodywash_progprob_effect;
            bodywash_progprob_effect = mean+pre_transformed_bodywash_progprob_effect;



            return bodywash_progprob_effect;


        }






//     ////effect on IQ due to secondary ISOLATION
//
//     double effect_of_secISO(gsl_rng * rng){
//     //effect of secondary isolation in reducing transmission from a colonized/infected isolated patient


//         double secISO_effect;
// //
//           secISO_effect= SECISOEFFECT;
//         return secISO_effect;
//     }


        double get_effect_of_secISO(gsl_rng * rng){

            double pre_transformed_secISO_effect;
            double sigma = effect_of_secISO_SD; //0.622;
            double mean = effect_of_secISO;//  0.365;
            pre_transformed_secISO_effect = gsl_ran_gaussian (rng, sigma);

            double secISO_effect;
            secISO_effect = mean+pre_transformed_secISO_effect;
            if(secISO_effect<0) secISO_effect=0;


            return secISO_effect;

        }


       int get_isocap(){

               return (ISOCAP)  ;
           }


///////////////////ECONOMICS DISTRIBUTIONS/////////////////////////

/*
    double get_bd_cost(gsl_rng * rng){
        double bd_cost;
        double alpha = 13.36;//for patients with 2 organs supported
        double beta = 100.12;//for patients with 2 organs supported
        bd_cost = gsl_ran_gamma (rng, alpha, beta);
        //cout<<bd_cost<<" bd_cost "<<"\n";
        return bd_cost;

    }

//mean for this is 1353

*/

///////////////////SENSITIVITY AND SPECIFICITY DISTRIBUTIONS/////////////////////////


/////conventional culture
        double get_CC_sensitivity(gsl_rng * rng){

            double pre_transformed_CC_sensitivity;
            double sigma = CC_SENSITIVITY_SD; // 0.194 by default;
            double mean = CC_SENSITIVITY_MEAN; // 0.6815 by default;
     pre_transformed_CC_sensitivity = gsl_ran_gaussian (rng, sigma);

            double CC_sens;
            CC_sens = mean+pre_transformed_CC_sensitivity;
            return CC_sens;
            cout<<CC_sens <<" CC_sens "<<"\n";

        }

        double get_CC_specificity(gsl_rng * rng){

            double pre_transformed_CC_specificity;
            double sigma =CC_SPECIFICITY_SD ;//  0.063 by default;
            double mean = CC_SPECIFICITY_MEAN ; //  0.882 by default;
            pre_transformed_CC_specificity= gsl_ran_gaussian (rng, sigma);
            double CC_spec;
            CC_spec = mean+pre_transformed_CC_specificity;
            return CC_spec;
        }

/////chromogenic agar (48h)
        double get_CA_sensitivity(gsl_rng * rng){

            double pre_transformed_CA_sensitivity;
            double sigma = CA_SENSITIVITY_SD; // defaults to 0.0427;
            double mean = CA_SENSITIVITY_MEAN ; //defaults to 0.8255;
            pre_transformed_CA_sensitivity = gsl_ran_gaussian (rng, sigma);

            double CA_sens;
            CA_sens = mean+pre_transformed_CA_sensitivity;

            return CA_sens;


        }



        double get_CA_specificity(gsl_rng * rng){

            double pre_transformed_CA_specificity;
            double sigma = CA_SPECIFICITY_SD;//  defaults to 0.1772;
            double mean =CA_SPECIFICITY_MEAN; //defaults to 0.8305;
            pre_transformed_CA_specificity= gsl_ran_gaussian (rng, sigma);

            double CA_spec;
            CA_spec = mean+pre_transformed_CA_specificity;

            return CA_spec;

        }

/////chromogenic agar early (24h)
        double get_CA_early_sensitivity(gsl_rng * rng){

            double pre_transformed_CA_early_sensitivity;
            double sigma = CA_EARLY_SENSITIVITY_SD;// defaults to 0.1249;
            double mean = CA_EARLY_SENSITIVITY_MEAN; //defaults to 0.6217;
            pre_transformed_CA_early_sensitivity = gsl_ran_gaussian (rng, sigma);

            double CA_early_sens;
            CA_early_sens = mean+pre_transformed_CA_early_sensitivity;

            return CA_early_sens;


        }



        double get_CA_early_specificity(gsl_rng * rng){

            double pre_transformed_CA_early_specificity;
            double sigma =  CA_EARLY_SPECIFICITY_SD; //defaults to 0.0417;
            double mean =CA_EARLY_SPECIFICITY_MEAN ; //defaults to 0.9713;
            pre_transformed_CA_early_specificity= gsl_ran_gaussian (rng, sigma);

            double CA_early_spec;
            CA_early_spec = mean+pre_transformed_CA_early_specificity;

            return CA_early_spec;


        }

/////PCR
        double get_PCR_sensitivity(gsl_rng * rng){

   //            double pre_transformed_PCR_sensitivity;
            double sigma = PCR_SENSITIVITY_SD; // 0.0510;
            double mean = PCR_SENSITIVITY_MEAN; //defaults to  0.8840;
            double pre_transformed_PCR_sensitivity = gsl_ran_gaussian (rng, sigma);

            double PCR_sens;
            PCR_sens = mean+pre_transformed_PCR_sensitivity;

            return PCR_sens;
        }


        double get_PCR_specificity(gsl_rng * rng){

            double pre_transformed_PCR_specificity;
            double sigma =  PCR_SPECIFICITY_SD; // defaults to  0.0474;
            double mean = PCR_SPECIFICITY_MEAN; //defaults to 0.8380;
            pre_transformed_PCR_specificity= gsl_ran_gaussian(rng, sigma);

            double PCR_spec;
            PCR_spec = mean+pre_transformed_PCR_specificity;

            return PCR_spec;
        }












