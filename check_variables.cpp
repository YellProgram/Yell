// A small program which just checks that variable definitions are ok
// it represents the part of yell input which Thomas now generates

//also get a handle to calculating the line number where error happens

//2:58
// checked how to read input arguments
// trying now to extract two simple functions - read file and check that file exists
// for this i am trying to isolate h5_basic_io.h
// failed
// isolated basic_io
// now can read file in. great.
//4:00
//4:16 медитирую над тем как это все устроить. Подключать всю грамматику из инпут файл парсера - лень. копипастить его нечестно
#include "basic_io.h"
#include "FormulaParser.h"
#include "precompiled_header.h"

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

using qi::lit;
using qi::char_;
using qi::eol;
using std::cout;
using std::string;

int main (int argc, char * const argv[]) {
    if (argc <= 1)
    {
      cout << "Usage: " << argv[0] << " <Filename>" << endl;
      exit(1);
    }
  
  string filename = string(argv[1]);
  
  if(!file_exists(filename))
  {
    cout << "Error: file '" << filename << "' does not exist";
    exit(1);
  }
  
  string file_content = read_input_file(filename); // Load the whole thing
  
  FormulaParser formula;
  
  Iterator start = file_content.begin();
  Iterator end = file_content.end();
  bool r;
  try{
    r = qi::phrase_parse(start,end,
                         qi::lexeme[formula], //qi::lexeme[(formula.assignment >> ';')]
                         ascii::space) ;//| lit("#") >> *(char_ - eol) >> eol
  }catch(const qi::expectation_failure<Iterator>& e)
  {
    string rest(e.first,e.last);
    cout << "ERROR: Parsing failed with exception\n" << rest;
    exit(1);
  }
  
  if(!r || start!=end)
  {
    string rest(start,end);
    cout << "ERROR: Parsing failed, the rest of the file is:\n" << rest;
    exit(1);
  }

  
  
  /*  skipper = boost::spirit::ascii::space | comment;
   comment = lit("#") >> *(char_ - eol) >> eol;
   
   refinable_parameters %=
   lexeme [ // for some unknown reasons skip(skipper_no_assignement) directive does not work. The program crashes with segfault. Thus I add skipper parser manually.
   lit("RefinableVariables")               > *skipper_no_assignement
   > "["                                   > *skipper_no_assignement
   >> *(valid_identifier                   > *skipper_no_assignement
   > "="                               > *skipper_no_assignement
   > double_                           > *skipper_no_assignement
   > -lit(';')                         > *skipper_no_assignement
   )                                   > *skipper_no_assignement
   > "]"
   ]
   */
  
  //skipper_no_assignement = boost::spirit::ascii::space | (lit("#") >> *(char_ - eol) >> eol);
  //skipper = omit[(formula.assignment >> ';')] | skipper_no_assignement;
  
 // comment = ;
  
  //Try to parse
  // check it is all right
  
  // Как насчет такого. Мы вручную добавим грамматику для скейла. перенесем set_refinable_parameters в грамматику калькулятора
  // У него как раз там будет логика чтобы инициализировать переменные.
  
  // выделим идентификатор и коммент в отдельный файл
  //
}


;

