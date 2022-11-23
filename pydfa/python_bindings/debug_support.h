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

#include <fstream>
#include <memory>
#include <boost/thread/shared_mutex.hpp>
#include <boost/thread.hpp>

#ifndef DEBUG_SUPPORT_H
#define DEBUG_SUPPORT_H

namespace math_model
{

namespace debug
{

extern std::ofstream debug_stream;

#ifdef DEBUG
extern boost::mutex debug_output_mutex;
#endif

#ifdef DEBUG
extern unsigned int debug_level;
#endif

void OpenDebugStream ();

void CloseDebugStream ();

std::ofstream & DebugLog ( unsigned int level );

template<class TYPE>
std::ofstream & DebugLog ( unsigned int level , TYPE message )
{
#ifdef DEBUG
	if ( ( debug_level >= level ) && ( debug_stream.is_open() ) )
	{
		debug_output_mutex.lock();
		debug_stream << "Level " << level << " : " << message << std::endl;
		debug_output_mutex.unlock();
	}
#endif
	return debug_stream;
}

template<class TYPE,class TYPE2>
std::ofstream & DebugLog ( unsigned int level , TYPE message , TYPE2 message2 )
{
#ifdef DEBUG
	if ( ( debug_level >= level ) && ( debug_stream.is_open() ) )
	{
		debug_output_mutex.lock();
		debug_stream << "Level " << level << " : " << message << " " << message2 << std::endl;
		debug_output_mutex.unlock();
	}
#endif
	return debug_stream;
}

template<class TYPE,class TYPE2, class TYPE3>
std::ofstream & DebugLog ( unsigned int level , TYPE message , TYPE2 message2 , TYPE3 message3)
{
#ifdef DEBUG
	if ( ( debug_level >= level ) && ( debug_stream.is_open() ) )
	{
		debug_output_mutex.lock();
		debug_stream << "Level " << level << " : " << message << " " << message2 << " " << message3 << std::endl;
		debug_output_mutex.unlock();
	}
#endif
	return debug_stream;
}

template<class TYPE,class TYPE2, class TYPE3, class TYPE4>
std::ofstream & DebugLog ( unsigned int level , TYPE message , TYPE2 message2 , TYPE3 message3 , TYPE4 message4)
{
#ifdef DEBUG
	if ( ( debug_level >= level ) && ( debug_stream.is_open() ) )
	{
		debug_output_mutex.lock();
		debug_stream << "Level " << level << " : " << message << " " << message2 << " " << message3 << " " << message4 << std::endl;
		debug_output_mutex.unlock();
	}
#endif
	return debug_stream;
}

template<class TYPE,class TYPE2, class TYPE3, class TYPE4, class TYPE5>
std::ofstream & DebugLog ( unsigned int level , TYPE message , TYPE2 message2 , TYPE3 message3 , TYPE4 message4 , TYPE5 message5)
{
#ifdef DEBUG
	if ( ( debug_level >= level ) && ( debug_stream.is_open() ) )
	{
		debug_output_mutex.lock();
		debug_stream << "Level " << level << " : " << message << " " << message2 << " " << message3 << " " << message4 << " " << message5 << std::endl;
		debug_output_mutex.unlock();
	}
#endif
	return debug_stream;
}

template<class TYPE,class TYPE2, class TYPE3, class TYPE4, class TYPE5, class TYPE6>
std::ofstream & DebugLog ( unsigned int level , TYPE message , TYPE2 message2 , TYPE3 message3 , TYPE4 message4 , TYPE5 message5 , TYPE6 message6)
{
#ifdef DEBUG
	if ( ( debug_level >= level ) && ( debug_stream.is_open() ) )
	{
		debug_output_mutex.lock();
		debug_stream << "Level " << level << " : " << message << " " << message2 << " " << message3 << " " << message4 << " " << message5 << " " << message6 << std::endl;
		debug_output_mutex.unlock();
	}
#endif
	return debug_stream;
}

}

}

#endif
