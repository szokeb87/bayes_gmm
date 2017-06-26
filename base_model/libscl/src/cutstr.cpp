/* ----------------------------------------------------------------------------

Copyright (C) 2004, 2006, 2008.

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

Function      cutstr - Cuts a string on white space or a delimiter.
                   

Syntax        #include "sclfuncs.h"
              std::string cutstr(char* str);
              std::string cutstr(std::string& str);
              std::vector<string> cutstr(const char* str, char delim);
              std::vector<string> cutstr(const std::string& str, char delim);

Prototype in  sclfuncs.h

Description   When called with one argument, str is either a null terminated 
              C string containing words separated by white space or a C++ 
              string containing words separated by white space.  On return,
              str contains the input string with the first word removed. 
            
              When called with two arguments, the first argument is either a 
              null terminated C string string containing words separated by 
              a delimiter or a C++ string containing words separated by a 
              delimiter.  The second argument is the delimiter.  There are 
              four cases.  
              1. The delimiter is white space other than a tab.
              2. The delimiter is a tab.
              3. The delimiter is a comma.
              4. The delimiter is other than the above.

Return value  When called with one argument, either the first word of str is 
              returned or a null string if str contains no words.  White
              space is stripped from the returned word.

              When called with two arguments, for each case as follows:

              1. A std::vector<string> containing the delimited words with
              white space stripped is returned.  If str is null or contains 
              only white space, the returned vector will be empty, i.e. will 
              have size zero.

              2. A std::vector<string> containing the delimited words is
              returned.  The tab is stripped.  White space is not stripped
              from the words.

              3. A std::vector<string> containing the delimited words is
              returned.  The comma is stripped.  White space is not 
              stripped from the words.  Quotes within the string are
              interpreted according to Excel CSV syntax.

              4. A std::vector<string> containing the delimited words is
              returned.  The delimiter is stripped.  White space is not 
              stripped from the words.

Remark        In case 3, if the the Excel syntax rules for quotes are 
              violated, the function will terminate with an error message
              and display the offending string.

Functions     Library: isspace
called        libscl: error

------------------------------------------------------------------------------*/

#include "sclfuncs.h"

std::string scl::cutstr(char* str)
{
  std::string remainder(str);
  std::string first_word;
  first_word = scl::cutstr(remainder);
  char* s = str;
  std::string::size_type n = remainder.size();
  for(std::string::size_type i=0; i<n; ++i) *s++ = remainder[i];
  *s = '\0';
  return first_word;
}

std::string scl::cutstr(std::string& str)
{
  const std::string null = "";
  std::string::size_type n = str.size();
  std::string::size_type i = 0;

  while (i<n && isspace(str[i])) ++i;       //find beginning first word
  if (i == n) {str = null; return null;}
  std::string::size_type begin_word = i;

  while (i<n && !isspace(str[i])) ++i;      //find one past end first word

  std::string first_word = str.substr(begin_word, i-begin_word);

  if (i == n) {str = null; return first_word;}

  str = str.substr(i,n-i);
  return first_word; 
}

std::vector<std::string> scl::cutstr(const char* str, char delim)
{
  return cutstr(std::string(str),delim);
}

std::vector<std::string> scl::cutstr(const std::string& str, char delim)
{
  const std::string null = "";
  std::vector<std::string> r;

  if (isspace(delim) && delim != '\t') { // Case 1.
    std::string s = str;
    std::string word = cutstr(s);
    if (word == null) return r;
    while (word != null) {
      r.push_back(word);
      word = cutstr(s);
    }
    return r;
  }

  std::string::size_type n = str.size();
  std::string::size_type i = 0;

  if (delim != ',') { // Case 2 or Case 4.
    while (i <= n) {
      if (i == n) {
        r.push_back(null); 
        return r;
      } 
      if (str[i] == delim) {
        r.push_back(null); 
        ++i;
      }
      else {
        std::string::size_type begin_word = i;
        while (i<n && str[i]!=delim) ++i;  // i equals n or str[i] is delim
        r.push_back(str.substr(begin_word, i-begin_word));
        if (i == n) {
          return r;
        } 
        else {
          ++i;
        }
      }
    }
    return r;
  }

  if (delim == ',') { // Case 3

    std::string s;
    s.reserve(n);
    std::string c;
    c.reserve(n);

    while (i < n) {
      if (str[i] == '"') {
        bool balanced_quotes = false;
        ++i;
        while (i < n) {
          if (i + 1 < n && str[i] == '"' && str[i+1] == '"') {
            s.push_back('"');
            c.push_back('"');
            i += 2;
          } 
          else if (str[i] == '"') {
            balanced_quotes = true;
            ++i;
            break;
          } 
          else if (str[i] == ',') {
            s.push_back('X');
            c.push_back(',');
            ++i;
          }
          else {
            s.push_back(str[i]);
            c.push_back(str[i]);
            ++i;
          }
        }
        if (i < n && str[i] != ',') balanced_quotes = false; 
        if ( !balanced_quotes ) {
          scl::error("Error, cutstr, unbalanced quotes in this line\n"+str);
        }
      }
      else {
        s.push_back(str[i]);
        c.push_back(str[i]);
        ++i;
      }
    }

    n = s.size();
    i = 0;

    while (i <= n) {
      if (i == n) {
        r.push_back(null); 
        return r;
      } 
      if (s[i] == delim) {
        r.push_back(null); 
        ++i;
      }
      else {
        std::string::size_type begin_word = i;
        while (i<n && s[i]!=delim) ++i;  // i equals n or s[i] is delim
        r.push_back(c.substr(begin_word, i-begin_word));
        if (i == n) {
          return r;
        } 
        else {
          ++i;
        }
      }
    }
    return r;
  }

  scl::error("Error, cutstr, this should never happen");
  return r;
}

