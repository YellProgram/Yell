//
// Created by Arkadiy on 19/02/2016.
//

#ifndef YELL_EXCEPTIONS_H
#define YELL_EXCEPTIONS_H

class TerminateProgram : std::exception {
public:
    TerminateProgram() {}
    ~TerminateProgram() throw() {}
};


#endif //YELL_EXCEPTIONS_H

