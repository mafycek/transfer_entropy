/***************************************************************************
 *   Copyright (C) 2005 by Hynek Lavička                                      *
 *   h.lavicka@email.cz                                                               *
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

#include <string.h>
#include <ostream>
#include <fstream>
#include <boost/thread/mutex.hpp>
#include "auxiliary_functions.h"

namespace math_model
{

namespace debug
{

#ifdef DEBUG
std::string debug_filename("debug.log");
#endif

std::ofstream debug_stream;

#ifdef DEBUG
boost::mutex debug_output_mutex;
#endif

#ifdef DEBUG
unsigned int debug_level = 8;
#endif

std::ofstream & DebugLog ( unsigned int level )
{
#ifdef DEBUG
	if ( ( debug_level >= level ) && ( debug_stream.is_open() ) )
	{
		debug_output_mutex.lock();
		debug_stream << "Level " << level << " : " ;
		debug_output_mutex.unlock();
	}
#endif
	return debug_stream;
}


void OpenDebugStream ()
{
#ifdef DEBUG
	std::string debug_full_filename;
	auxiliary_functions:: MakeFileName ( debug_full_filename , debug_filename );
	debug_stream.open(debug_full_filename.data());
	std::locale::global(std::locale(""));
	debug_stream.imbue(std::locale());
#endif
}

void CloseDebugStream ()
{
#ifdef DEBUG
	debug_stream.close();
#endif
}

}

}
