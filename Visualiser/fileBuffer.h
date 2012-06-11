#pragma once

#include <string>
#include <fstream>      //allows file io
class fileBuffer {
private:
bool firstLine;
std::ifstream logFile;

public:
double length;
int endFile;
  void openFile()
  {
    logFile.open("locLog.dat");
  }
  void bufferFile(int bufferSize)
  {
    events.clear();
    endFile = 100;
    for(int i = 0; i < bufferSize; ++i)
    {
      if(!logFile.good())
      {
        endFile = i;
        break;
      }

      std::string line;
      getline(logFile, line); //get next line of log
      processLine(line);

    }

  }
  void processLine(std::string line)
  {
    int counter = 0;
    double x = 0, y = 0, z = 0;
    double time = 0;
    int strpos = 1;

  std::vector<Particle> p;
    if(firstLine)
      ++counter;
    for(int i = 0; i < line.size(); ++i)
    {
      if(line[i] == '\t')
      {
        switch(counter)
        {
          case 0: //then system time
            time = atof(line.substr(0, i).c_str());
            strpos = i + 1;
            break;
          case 1: //then x
            x = atof(line.substr(strpos,i - strpos).c_str());
            strpos = i + 1;
            break;
          case 2: //then y
            y = atof(line.substr(strpos,i - strpos).c_str());
            strpos = i + 1;
            break;
          case 3: //then z
            z = atof(line.substr(strpos,i - strpos).c_str());
            strpos = i + 1;
            counter = 0; //get ready for x from next particle
            if(!firstLine)
              p.push_back(Particle(x,y,z)); //add to particles
            else
              length = y;
            break;
          default:
            exit(1);
            std::cout << "error in reading file" << std::endl;
        }
        ++counter;
      }
    }
    if(firstLine)
      firstLine = false;
    else
      events.push_back(Event(time, p));
  }

};


