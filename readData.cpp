#include "main.h"

#define N_GEOM 4
#define N_PARAM 2

int readData(std::string fileName, Truss &truss, Parameters &parameters){
    // Get the file
    std::ifstream inFile(fileName);
    std::string buf;
    std::getline(inFile, buf);
    // Check the file type
    std::string fileType = "#PARAMETERS";
    assert(buf == fileType);
    // Loop on the file sections
    while(true){
        std::getline(inFile, buf);
        buf.erase(std::remove(buf.begin(), buf.end(), ' '),buf.end());
        switch(buf[0]){
            case '%' : // Comment identifier
            continue;
            break;
            case '#' : // Starting parameter list (need to be ordered)
            buf.erase (0,1);
            if(buf=="geom"){
                int cnt=0;
                char value[40];
                while(cnt < N_GEOM){
                    std::getline(inFile, buf);
                    if(1==sscanf(buf.c_str(),"%*[^=]=%s", value)){
                        switch(cnt){ // Parameter filling
                            case 0:
                            truss.A = atof(value);
                            break;
                            case 1:
                            truss.a = atof(value);
                            break;
                            case 2:
                            truss.b = atof(value);
                            truss.l0 = sqrt(truss.a*truss.a + truss.b*truss.b);
                            break;
                            case 3:
                            truss.E = atof(value);
                            break;
                        }
                        ++cnt;
                    }
                    else{
                        std::cout << "Not enough arguments or white line" << std::endl;
                        return 1;
                    }
                }
            }
            else if(buf=="param"){
                int cnt=0;
                char value[40];
                while(cnt < N_PARAM){
                    std::getline(inFile, buf);
                    if(1==sscanf(buf.c_str(),"%*[^=]=%s", value)){
                        switch(cnt){ // Parameter filling
                            case 0:
                            parameters.resolution = atoi(value);
                            break;
                            case 1:
                            parameters.tFinal = atof(value);
                            break;
                        }
                        ++cnt;
                    }
                    else{
                        std::cout << "Not enough arguments or white line" << std::endl;
                        return 1;
                    }
                }
            }
            else if(buf=="end"){break;}
            else if(buf=="END"){return 0;}
            else{
                std::cout <<"Unknown '"<<buf<<"' identifier."<< "\n";
                return 1;
            }
            break;
            default :
            continue;
        }
    }
}
