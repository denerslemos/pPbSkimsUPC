#include "call_libraries.h"  // call libraries from ROOT and C++

// print UIC Jets welcome message
/*
Arguments
printout: true print UIC Jet logo, false do not print out
*/
void printwelcome(bool printout){
	if(!printout) return;
	cout << endl;
	cout << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "         << endl;
	cout << "+                                                              +"         << endl;
	cout << "+                                      ______        ______    +"         << endl;
	cout << "+   UUU    UUU  IIIIIIII      CCCCCC  |   % /   /\\   \\ %   |   +"         << endl;
	cout << "+   UUU    UUU    IIII      CCCC      |  % /   /  \\   \\ %  |   +"         << endl;
	cout << "+   UUU    UUU    IIII    CCC         | % /   / %  \\   \\ % |   +"         << endl;
	cout << "+   UUU    UUU    IIII    CCC         |% /   / Jets \\   \\ %|   +"         << endl;
	cout << "+   UUU    UUU    IIII      CCCC      | /   /   %    \\   \\ |   +"         << endl;
	cout << "+    UUUUUUUU   IIIIIIII      CCCCCC  |/   /__________\\   \\|   +"         << endl;
	cout << "+                                                              +"         << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"         << endl;
	cout << "+                                                              +"         << endl;
	cout << "+    Welcome to UIC-HENP pPb skim code for Data and MC         +"         << endl;
	cout << "+            From Dener Lemos: dener.lemos@cern.ch             +"         << endl;
	cout << "+                                                              +"         << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"         << endl;
}

// print start day and time message
void print_start(){
	cout << endl;
    time_t init = time(0);
    char* init_time = ctime(&init); // convert now to string form
    cout << "Starting at : " << init_time << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
	cout << endl;
}

// print stop day and time message
void print_stop(){
    time_t end = time(0);
    char* end_time = ctime(&end); // convert now to string form
   	cout << endl;
	cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   	cout << endl;
	cout << "Stopping at : " << end_time << endl;
	cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << endl;
}
