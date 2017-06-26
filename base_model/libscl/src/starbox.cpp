/*-----------------------------------------------------------------------------

Copyright (C) 1990, 1993, 1994, 2002, 2003, 2006, 2012.

A. Ronald Gallant
Post Office Box 659
Chapel Hill NC 27514-0659
USA   

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

-------------------------------------------------------------------------------

Class         starbox  - Puts title in a centered box of stars with usage 
                         in the style of a C++ manipulator; e.g.,
                         cout << starbox("/title//"); or
                         cout << starbox("\n title \n\n",'\n');

Syntax        #include "sclfuncs.h"
              
              class starbox {
              private:
                char boxed_title[MAX_STARBOX_SIZE];
              public:
                starbox();
                starbox(const starbox& sb);
                starbox(const std::string& title, char crlf = '/');
                starbox(const char *title, char crlf = '/');
                starbox& operator=(const starbox& sb);
                string get_boxed_title() const;
                string operator()() const;      // short for get_boxed_title()
                friend std::ostream& operator<<(ostream& os, const starbox& sb);
              };

Description   title is a pointer to a string literal, a pointer to an array 
              of char, or a std::string.  The string pointed to or referenced 
	      by title is made up of concatonated lines.  A line is made up 
	      of 0-68 characters followed by a slash (default). Examples are:
              cout << starbox("//line 1/line 2/line 3///");
              cout << starbox("\n\nline 1\nline 2\nline 3\n\n\n",'\n');

Functions     Library: strcat, strchr, strlen
called        Libcpp:   (none)

-----------------------------------------------------------------------------*/

#include "sclfuncs.h"

using std::strcat;
using std::strchr;
using std::strlen;

namespace scl {

  void starborder(const char *title, size_t sb_len, char* sb_str, char crlf);
  
  starbox::starbox()
  { 
    const char* title = "/no title//";
    starborder(title, MAX_STARBOX_SIZE, boxed_title, '/');
  }

  starbox::starbox(const char *title, char crlf)
  {
    starborder(title, MAX_STARBOX_SIZE, boxed_title, crlf);
  }

  starbox::starbox(const std::string& title, char crlf)
  {
    starborder(title.c_str(), MAX_STARBOX_SIZE, boxed_title, crlf);
  }

  starbox::starbox(const starbox& sb)
  {
    strcpy(boxed_title, sb.boxed_title);
  }
  
  starbox& starbox::operator=(const starbox& sb)
  {
    if (this != &sb) {
      strcpy(boxed_title, sb.boxed_title);
    }
    return *this;
  }

  std::string starbox::get_boxed_title() const
  {
    return std::string(boxed_title);
  }
  
  std::string starbox::operator()() const
  {
    return std::string(boxed_title);
  }
  
  std::ostream& operator<<(std::ostream& os, const starbox& sb)
  {
    os << sb.boxed_title;
    return os;
  }
  
  void starborder(const char *title, size_t sb_len, char* sb_str, char crlf)
  {
    const char* begin;
    const char* end;
    char        border[81], tab[81], work[256], middle[81];
    char*       line;
    char*       t;
    int         i, pad, mpad, linesize;
    int         length;
    size_t      cum_length;
  
    const char* err_msg = "\n"
         "\t**********************************************\n"
         "\t*      Error, starbox, title too long        *\n" 
         "\t**********************************************\n";

    cum_length = 0;
    *sb_str = '\0';
  
    for (i=0; i<70; i++)  border[i] = '*';
    border[70] = '\n'; border[71] = '\0';
  
    ((LINESIZE<72)||(LINESIZE>133)) ? (linesize=133) : (linesize=LINESIZE);
    
    pad = (linesize-72)/2+1;
    for (i=0; i<pad; i++)  tab[i] = ' ';
    tab[pad] = '\0';
  
    work[0]='\n'; work[1]='\n'; work[2]='\0';
    line = strcat(work,tab);
    line = strcat(work,border);
  
    cum_length += strlen(line);
    if (cum_length < sb_len) {
      strcat(sb_str,line);
    }
    else {
      strcpy(sb_str,err_msg);
      return;
    } 

    begin = title;
    while ( (end=strchr(begin,crlf)) != 0 ) {

      length = end - begin;
      if (length > 68) {
        strcpy(sb_str,err_msg);
        return;
      }

      middle[0] = '*';
      mpad = (68-length)/2+1;
      for (i=1; i<=mpad; i++) middle[i] = ' ';
      t = &middle[mpad];
      while (begin < end) *t++ = *begin++;
      begin++;
      for (i=mpad+length; i<=68; i++) middle[i] = ' ';
      middle[69]='*'; middle[70]='\n'; middle[71]='\0';
      work[0]='\0';
      line = strcat(work,tab);
      line = strcat(work,middle);
  
      if (cum_length < sb_len) {
        strcat(sb_str,line);
      }
      else {
        strcpy(sb_str,err_msg);
        return;
      }
    }
  
    work[0]='\0';
    line = strcat(work,tab);
    line = strcat(work,border);
  
    if (cum_length < sb_len) {
      strcat(sb_str,line);
    }
    else {
      strcpy(sb_str,err_msg);
      return;
    }
  
    return;
  
  }
  
} 
