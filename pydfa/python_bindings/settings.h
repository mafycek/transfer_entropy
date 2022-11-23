/*! \file settings.h
\brief A general macro and variable declaring file

You can find here all general setting of the program
*/
/***************************************************************************
 *   Copyright (C) 2003 by Hynek Lavička                                   *
 *   h.lavicka@email.cz                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef MATH_SETTINGS_H
#define MATH_SETTINGS_H
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <climits>

#ifndef ULLONG_MAX
#define ULLONG_MAX 18446744073709551615ULL
#endif

/*!
 \brief type which will be used as identificator of nodes, agents and other objects.
 */
#define TYPE_OF_NODE unsigned long

/*!
 \brief Define how many agents can be bound to some node with a bind.
 */
#define NODE_WITH_BIND_MAX_NUMBER 8

/*!
 \brief Define how exponent for utility function.
 */
#define EXPONENT_FOR_UTILITY_FUNCTION 2

/*!
 \brief Define that is consisted of name of base output directory.
Deprecated symbol.
 */
#define BASE_OUTPUT_DIRECTORY std::string( std::string ( getenv ( "HOME" ) ) + "/work/computations/")
//#define BASE_OUTPUT_DIRECTORY "/mnt/hda9/guest/results/"

/*!
 \brief Define that is consisted of name file where debug output will be stored.
 */
#define DEBUG_OUTPUT "debug.log"

/*!
 \brief Define that is consisted of name of general simulation information log.
 */
#define GENERAL_INFO_OUTPUT std::string("general_simulation_information.log")
