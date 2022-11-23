/*! \file auxiliary_functions.h
\brief Auxilliary functions declaration file.

Auxilliary functions are function that are only helping function.
\ingroup zarja
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

#ifndef AUXILLIARY_FUNCTIONS_H
#define AUXILLIARY_FUNCTIONS_H

#include <string>
#include <list>

#include "settings.h"

/**
\brief Namespace encapulates all auxilliary functions.
\ingroup zarja
*/
namespace auxiliary_functions
{

/// A function which takes no arguments and it gets an integer
/**
The function returns output directory
\return output directory
*/
std::string & GetOutputDir ( );


/// A function which takes one argument and it gets an integer
/**
The function sets up output dir
\param dir New output directory
\return 0 if all OK otherwise non-zero value
*/
int SetOutputDir ( const std::string & dir );


/// A function which takes no arguments
/**
The function sets up output dir to relative value.

\return 0 if all OK otherwise non-zero value
*/
int SetOutputDir ( );

/// A function which takes no arguments and it gets an integer
/**
The function returns relative output directory
\return relative output directory
*/
std::string & GetRelativeOutputDir ( );


/// A function which takes one argument and it gets a pointer to char
/**
The function sets relative output directory
\param dir String where is new directory located
\return If all OK than 0, otherwise non-zero value.
*/
int SetRelativeOutputDir ( const char * dir );


/// A function which takes no arguments
/**
The function gets OutputBaseDir variable.

\return 0 if all OK otherwise non-zero value
*/
std::string & GetOutputBaseDir ( );


/// A function which takes no arguments
/**
The function sets up OutputBaseDir default value.

\return 0 if all OK otherwise non-zero value
*/
int SetOutputBaseDir ( );


/// A function which takes 4 arguments and it gets an integer
/**
The function gets filename which is consisted of general_dir, specific_dir and file
\param output String where will be the result stored
\param general_dir A general dir of the project
\param specific_dir A specific dir of the project
\param file Name of file
\return If all OK than 0, otherwise non-zero value.
*/
int MakeFileName ( std::string & output , const std::string & general_dir , const std::string & specific_dir , const std::string & file );


/// A function which takes 2 arguments and it gets an integer
/**
The function gets filename which is consisted of file
\param output String where will be the result stored
\param file Name of file
\return If all OK than 0, otherwise non-zero value.
*/
int MakeFileName ( std::string & output , const std::string & file );


/// A function which takes 3 arguments and it gets an integer
/**
The function gets filename which is consisted of file
\param output String where will be the result stored
\param file Name of file
\param suffix Suffix of the file
\return If all OK than 0, otherwise non-zero value.
*/
int MakeFileName ( std::string & output , const std::string & file , const std::string & suffix );


/// A function which takes 3 arguments and it gets an integer
/**
The function gets filename which is consisted of file
\param output String where will be the result stored
\param file Name of file
\param number Number of the file
\param suffix Suffix of the file
\return If all OK than 0, otherwise non-zero value.
*/
int MakeFileName ( std::string & output , const std::string & file , int number , const std::string & suffix );


/// A function which takes 4 arguments and gets an integer.
/**
The function gets minimum and maximum of the data stored in file
\param input Name of the file which will be processed
\param minimum Place where to store minimum, nullptr drops the value.
\param maximum Place where to store maximum, nullptr drops the value.
\param rows Place where to store number of rows, nullptr drops the value.
\return If all OK than 0, otherwise non-zero value.
*/
int GetExtremesFileData1stRowInteger ( char * input , int * minimum , int * maximum , int * rows );


/// A function which takes 3 arguments a pointer to char and 2 doubles and 1 integer and gets an integer.
/**
The function gets minimum and maximum of the data stored in file
\param input Name of the file which will be processed
\param minimum Place where to store minimum, nullptr drops the value.
\param maximum Place where to store maximum, nullptr drops the value.
\param rows Place where to store number of rows, nullptr drops the value.
\return If all OK than 0, otherwise non-zero value.
*/
int GetExtremesFileData1stColumnDouble ( char * input , double * minimum , double * maximum , int * rows);


/// A function which takes 4 arguments and gets an integer.
/**
	The function gets minimum and maximum of the data stored in file
	\param input Name of the file which will be processed
	\param minimum Place where to store minimum, nullptr drops the value.
	\param maximum Place where to store maximum, nullptr drops the value.
	\param rows Place where to store number of rows, nullptr drops the value.
	\return If all OK than 0, otherwise non-zero value.
 */
int GetExtremesFileData1stColumnInteger ( char * input , int * minimum , int * maximum , int * rows );


/// A function which takes 3 arguments and gets an integer.
/**
	The function load an array with data from first column of file.
	\param input Name of the file which will be processed
	\param array Pointer where are stored data.
	\param elements Pointer where number of loaded data will be stored.
	\return If all OK than 0, otherwise non-zero value.
 */
int Load1stColumnDoubleToArray ( char * input , double * * array , int * elements );


/// A function which takes 3 arguments and gets an integer.
/**
	The function destroy array which is allocated using function Load1stColumnDoubleToArray.
	\param array Pointer to data.
	\return If all OK than 0, otherwise non-zero value.
 */
int DestroyDoubleArray ( double * array );


int compare_TYPE_OF_NODE ( const void * a , const void * b );

/// A function with 2 arguments that gets two double value
/**
	The function compare two double values and returns result.
	\param pointer1 The first double.
	\param pointer2 The second double.
	\return -1 when a > b, 0 when a = b and 1 when a < b.
 */
int compare_double ( const void * pointer1 , const void * pointer2 );

/// A function with 2 arguments that gets two integer value
/**
	The function compare two integer values and returns result.
	\param pointer1 The first integer.
	\param pointer2 The second integer.
	\return -1 when a > b, 0 when a = b and 1 when a < b.
 */
int compare_int ( const void * a , const void * b );

/// A function with 2 arguments that gets two complex<double> value
/**
	The function compare two complex<double> values by real part and returns result.
	\param pointer1 The first integer.
	\param pointer2 The second integer.
	\return -1 when a > b, 0 when a = b and 1 when a < b.
 */
int compare_complex_double ( const void * a , const void * b );

void RegisterToCollector ( void * memory );

void CleanUpCollector ( );

}

#endif
