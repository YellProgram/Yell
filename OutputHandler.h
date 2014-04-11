//
//  OutputHandler.h
//  diffuser_y
//
//  Created by Arkadiy Simonov on 2/26/13.
//
//

#ifndef __diffuser_y__OutputHandler__
#define __diffuser_y__OutputHandler__

#include <iostream>
#include <string>

using std::string;
using std::ostream;
using std::cout;

#define REPORT(when) if(report(when)) cout 

enum WhenToPrint { FIRST_RUN, AFTER_REFINEMENT, MAIN, ERROR };

/**
 *   Class to output the messages in a standardized way. For now just outputs a bool which tells whether to print or not the message. The class atm is used through macros
 */
class OutputHandler {
public:
  OutputHandler()
  {
    first_run_flag=false;
    last_run_flag=false;
    shut_up_flag=false;
  }
  
  /**
   Flag to output the information from the first model construction
   this includes echoing the information on say the number of pairs
   created.
   */
  bool first_run_flag;
  /**
   Flag to output the information about model after the refinement is finished.
   eg. print out values of variables
   */
  bool last_run_flag;
  
  /**
   Flag for test purposes. Suppresses any output
   */
  bool shut_up_flag;
  
  void last_run()  {
    last_run_flag = true;
  }
  
  /**
   Supresses any further output.
   */
  void shut_up()  {
    shut_up_flag=true;
  }
  
  void expect_first_calculation() {
    first_run_flag=true;
  }
  
  void calculation_is_finished()  {
    first_run_flag = false;
    last_run_flag = false;
  }
  
  /** The function that handles all the output from the program.
   
  The program execution is not extremely nice at the moment, since most of the logic is executed while the input file gets parsed (which happens for every calculation). This class helps to keep track on what kind of output should be printed when.
   */
  bool operator()(WhenToPrint when) {
    if (when==ERROR)
    {
      cout << "\n ERROR: ";
      return true;
    }
    
    if (shut_up_flag)
      return false;
    
    switch (when) {
      case FIRST_RUN:
        return first_run_flag;
        break;
      case AFTER_REFINEMENT:
        return last_run_flag;
        break;
      case MAIN:
        return true;
        break;
      case ERROR:
        return true;
        break;
    }
    return true;
  }
};

#endif /* defined(__diffuser_y__OutputHandler__) */
