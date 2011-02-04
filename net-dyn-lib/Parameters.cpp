#include "Parameters.h"
#include <fstream>
#include <list>
#include <algorithm>
#include <functional>
#include <libgen.h>

void Parameters::handleArgs(int argc, char **argv)
{ static string dirname_for_settings;
  dirname_for_settings = dirname(argv[0]);
  cout << "Executable directory is " << dirname_for_settings << endl;
  for ( ++argv; (*argv); ++argv )
  {
    string arg(*argv);
    if (arg=="-f")
    {
      //cout << "argument: -f " << argv[1] << endl;
      string filename(*++argv);
      parseSettingsFile(filename);
    }
    else if (arg[0] == '-' && arg[1] == '-')
    { // this currently can only be --param=val
      string line(arg,2,string::npos);
      string::size_type eq = line.find('=');
      if (eq != string::npos)
        line[eq] = ' ';
      //cout << line << endl;
      parseLine(line);
    }
    else
    {
      cout << "unknown argument: " << arg << endl;
    }
  }
}

void Parameters::parseSettingsFile(string filename)
{ static string dirname_for_settings;
  string pathname = dirname_for_settings + '/' + filename;
  cout << "parseSettingsFile(\"" << pathname << "\")\n";
  ifstream settings(pathname.c_str());
  if (settings.is_open())
    parseSettings(settings);
  else
  { cout << "couldn't open " << pathname << endl;
    exit(-1);
  }
}

// each line here must be
// "key val"
// with just one space between the two;
// or, "#comment"; or, "include filename";
// otherwise the line will be ignored or misinterpreted.
// leading and trailing spaces are ignored.
void Parameters::parseSettings(istream &file)
{
  string line;
  while (getline(file,line))
  {
    parseLine(line);
  }
}

string Parameters::cleanLine(string line)
{
  // get rid of anything starting with #
  string::size_type pound = line.find('#');
  if (pound != string::npos)
    line.erase(pound);

  // get rid of leading + trailing space
  string::iterator e = line.end();
  while (e != line.begin() && isspace(e[-1]))
    --e;
  string::iterator b = line.begin();
  while (b != e && isspace(*b))
    ++b;
  //line.erase(e,line.end());
  line.replace(line.begin(),line.end(),b,e);

  return line;
}

void Parameters::parseLine(string line)
{
  line = cleanLine(line);

  // check for include directive
  string directive("include ");
  if (line.substr(0,directive.size())==directive)
  {
    parseSettingsFile(line.substr(directive.size(),string::npos));
    return;
  }

  // what's left has to have a space as separator
  string::size_type space = line.find(' ');
  if (space == string::npos)
    return;
  string key(line,0,space);
  string val(line,space+1);
  cout << key << "=" << val << endl;
  set(key,val);
}

void Parameters::writeAllSettings(ostream&outs)
{
  list<string> keys;
  // get all keys from the dictionary
  transform(begin(),end(), inserter(keys,keys.begin()),
	    _Select1st< pair<string,string> >());
  // output key/value pairs using sorted order
  for (list<string>::iterator ki = keys.begin(); ki != keys.end(); ++ki)
    if (const string *getki = get(*ki))
      outs << *ki << ' ' << *getki << '\n';
}
