/*
 *  precompiled_header.h
 *  diffuser_y
 *
 *  Created by Arkadiy Simonov on 3/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_real.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/home/phoenix/bind/bind_member_function.hpp>
#include <boost/spirit/home/phoenix/operator/member.hpp>

#include <boost/spirit/home/phoenix/operator/self.hpp>
#include <boost/spirit/home/phoenix/object/construct.hpp>
#include <boost/spirit/home/phoenix/bind/bind_function.hpp>
#include <boost/spirit/home/phoenix/object/new.hpp>
#include <boost/variant.hpp>
#include "boost/tuple/tuple.hpp"

#include <vector>
#include <string>

namespace phoenix = boost::phoenix;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

typedef std::string::iterator iterator_ ;