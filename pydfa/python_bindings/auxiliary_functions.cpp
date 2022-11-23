/*! \file auxiliary_functions.cpp
\brief Auxilliary functions definition file.

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

#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <complex>
#include "auxiliary_functions.h"

namespace auxiliary_functions
{

/// Static string.
/**
	The variable stores relative name (with respect OutputDir) of directory which will be used for assembling output directory name.
 */
static std::string RelativeOutputDir;

/// Static string.
/**
The variable stores whole name of directory where all outputs will be produced.
*/
static std::string OutputDir;

/// Static string.
/**
The variable stores name of base name of directory where all output will be produced.
*/
static std::string OutputBaseDir;

std :: list<void *> zarja_garbage_collector;

std::string & GetRelativeOutputDir ( )
{
    return RelativeOutputDir;
}


int SetRelativeOutputDir ( const char * dir )
{
    RelativeOutputDir = std::string ( dir );

    return 0;
}


std::string & GetOutputDir ( )
{
    return OutputDir;
}


int SetOutputDir ( )
{
    SetOutputBaseDir ( );
    OutputDir = GetOutputBaseDir() + "/" + GetRelativeOutputDir ( ) + "/";

    return 0;
}


int SetOutputDir ( const std::string & dir )
{
    OutputDir = dir +GetRelativeOutputDir ( ) + "/";

    return 0;
}


std::string & GetOutputBaseDir ( )
{
    return OutputBaseDir;
}


int SetOutputBaseDir ( )
{
    char * home;
    home = getenv ("HOME");

    if ( home == nullptr )
    {
        std :: cout << "SetOutputBaseDir: Cannot get HOME enviroment variable." << std :: endl;

        return 1;
    }
    else
    {
        OutputBaseDir = std::string ( home ) + "/work/computations";

        return 0;
    }
}


int MakeFileName ( std::string & output , const std::string & general_dir , const std::string & specific_dir , const std::string & file )
{
    output = general_dir + specific_dir + "/" + file;

    return 0;
}


int MakeFileName ( std::string & output , const std::string & file )
{
    output = GetOutputDir ( ) + file;

    return 0;
}


int MakeFileName ( std::string & output , const std::string & file , const std::string & suffix )
{
    output = GetOutputDir ( ) + file + suffix;

    return 0;
}


int MakeFileName ( std::string & output , const std::string & file , int number , const std::string & suffix )
{
    char number_char [ 100 ];

    sprintf ( number_char , "%i" , number );

    output = GetOutputDir ( ) + "/" + file + "-" + number_char + suffix;

    return 0;
}


int GetExtremesFileData1stRowInteger ( char * input , int * minimum , int * maximum , int * rows)
{
    FILE * input_file;
// opening output stream
    if ( nullptr == ( input_file = fopen ( input , "r" ) ) )
    {
        return 1;
    }

// analyzing data
    long counter = 0;
    int min = INT_MAX , max = INT_MIN;
    for ( ; ! feof ( input_file ) ; )
    {
        int event;
        char character = 65;
        fscanf ( input_file , "%i" , & event );
// searching for new line
        for ( ; ! feof ( input_file ) || character != '\n' ; )
        {
            character = ( char ) fgetc ( input_file );
        }

        if ( ! feof ( input_file ) )
        {
            counter ++;
        }
// searching for new maximum and minimum
        if ( event > max )
        {
            max = event;
        }
        if ( event < min )
        {
            min = event;
        }
    }

// close stream
    fclose ( input_file );

// value transaction
    if ( minimum != nullptr )
    {
        * minimum = min;
    }

    if ( maximum != nullptr )
    {
        * maximum = max;
    }

    if ( rows == nullptr )
    {
        * rows = counter;
    }

    return 0;
}


int GetExtremesFileData1stColumnInteger ( char * input , int * minimum , int * maximum , int * rows )
{
    FILE * input_file;
// opening output stream
    if ( nullptr == ( input_file = fopen ( input , "r" ) ) )
    {
        return 1;
    }

// analyzing data
    long counter = 0;
    int min = INT_MAX , max = INT_MIN;
    for ( ; ! feof ( input_file ) ; )
    {
        int event;
        char character = 65;
        fscanf ( input_file , "%i" , & event );
// searching for new line
        for ( ; ! feof ( input_file ) && character != '\n' ; )
        {
            character = ( char ) fgetc ( input_file );
        }

        if ( ! feof ( input_file ) )
        {
            counter ++;
        }
// searching for new maximum and minimum
        if ( event > max )
        {
            max = event;
        }
        if ( event < min )
        {
            min = event;
        }
    }

// close stream
    fclose ( input_file );

// value transaction
    if ( minimum != nullptr )
    {
        * minimum = min;
    }

    if ( maximum != nullptr )
    {
        * maximum = max;
    }

    if ( rows != nullptr )
    {
        * rows = counter;
    }

    return 0;
}

int GetExtremesFileData1stColumnDouble ( char * input , double * minimum , double * maximum , int * rows )
{
    FILE * input_file;
// opening output stream
    if ( nullptr == ( input_file = fopen ( input , "r" ) ) )
    {
        return 1;
    }

// analyzing data
    long counter = 0;
    double min = 1.7e308 , max = -1.7e308;
    for ( ; ! feof ( input_file ) ; )
    {
        double event;
        char character = 65;
        fscanf ( input_file , "%lg" , & event );
        // searching for new line
        for ( ; ! feof ( input_file ) && character != '\n' ; )
        {
            character = ( char ) fgetc ( input_file );
        }

        if ( ! feof ( input_file ) )
        {
            counter ++;
        }
// searching for new maximum and minimum
        if ( event > max )
        {
            max = event;
        }
        if ( event < min )
        {
            min = event;
        }
    }

// close stream
    fclose ( input_file );

// value transaction
    if ( minimum != nullptr )
    {
        * minimum = min;
    }

    if ( maximum != nullptr )
    {
        * maximum = max;
    }

    if ( rows != nullptr )
    {
        * rows = counter;
    }

    return 0;
}


int Load1stColumnDoubleToArray ( char * input , double * * array , int * elements )
{
    int rows;
    auxiliary_functions :: GetExtremesFileData1stColumnDouble ( input , nullptr , nullptr , & rows);
    double * data_array = ( double * ) malloc ( sizeof ( double ) * ( rows + 1 ) );

    if ( data_array == nullptr )
    {
        return 1;
    }

    FILE * input_file;

    if ( nullptr == ( input_file = fopen ( input , "r" ) ) )
    {
        return 1;
    }

    int counter = 0;
    for ( ; ! feof ( input_file ) ; )
    {
        double event;
        char character = 65;
        fscanf ( input_file , "%lg" , & event );
        data_array [ counter ] = event;
// searching for new line
        for ( ; ! feof ( input_file ) && character != '\n' ; )
        {
            character = ( char ) fgetc ( input_file );
        }
        if ( ! feof ( input_file ) )
        {
            counter ++;
        }
    }

    * array = data_array;

    if ( elements != nullptr )
    {
        * elements = counter;
    }

    fclose ( input_file );
    return 0;
}


int DestroyDoubleArray ( double * array )
{
    if ( array != nullptr )
    {
        free ( array );
        return 0;
    }
    return 1;
}


int compare_TYPE_OF_NODE ( const void * a , const void * b )
{
    return ( ( TYPE_OF_NODE * ) a ) [ 0 ] > ( ( TYPE_OF_NODE * ) b ) [ 0 ];
}


int compare_double ( const void * a , const void * b )
{
    if ( ( * ( double * ) a ) > ( * ( double * ) b ) )
    {
        return -1;
    }
    if ( ( * ( double * ) a ) == ( * ( double * ) b ) )
    {
        return 0;
    }
    if ( ( * ( double * ) a ) < ( * ( double * ) b ) )
    {
        return 1;
    }
    return 0;
}


int compare_int ( const void * a , const void * b )
{
    return ( ( * ( int * ) a ) - ( * ( int * ) b ) );
}


int compare_complex_double ( const void * a , const void * b )
{
    if ( ( * ( std::complex<double> * ) a ) . real ( ) > ( * ( std::complex<double> * ) b ) . real () )
    {
        return -1;
    }
    if ( ( * ( std::complex<double> * ) a ) . real ( ) == ( * ( std::complex<double> * ) b ) . real ( ) )
    {
        return 0;
    }
    if ( ( * ( std::complex<double> * ) a ) . real ( ) < ( * ( std::complex<double> * ) b ) . real ( ) )
    {
        return 1;
    }
    return 0;
}


void RegisterToCollector ( void * memory )
{
    zarja_garbage_collector.push_back ( memory );
}


void CleanUpCollector ( )
{
    std :: list<void *>::iterator end_iter = zarja_garbage_collector.end( );
    int i = 0;

    for ( std :: list<void *>::iterator iter = zarja_garbage_collector.begin( ) ; iter != end_iter ; ++ iter )
    {
        i ++;
        void * memory_strip = ( * iter );
// 			delete memory_strip;
    }

    zarja_garbage_collector.clear ( );
}

}



