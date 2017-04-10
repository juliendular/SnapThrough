#include "main.h"

int writeData(std::vector<double> &x, std::vector<double> &y, std::string fileName){
    // Check the vector lengths
    if(x.size() != y.size()){
        std::cout << "Vectors must have the same length" << std::endl;
        return 1; // Error code
    }
    // Create the file (or erase an existing one)
    fileName = "results/" + fileName + ".txt";
    const char *fileNameChar = fileName.c_str();
    FILE *fp = fopen(fileNameChar, "w");
    if(fp != NULL){
        // Write the x values
        for(int i=0 ; i < x.size() ; i++)
            fprintf(fp, "%f ", x[i]);
        fprintf(fp, "\n");
        // Write the y values
        for(int i=0 ; i < y.size() ; i++)
            fprintf(fp, "%f ", y[i]);
        fclose(fp);
        return 0;
    }
    else{
        printf("Cannot create a .txt file\n");
        return 2; // Error code
    }
}

int writeData(std::vector<std::vector<double> > &x, std::vector<double> &y, std::string fileName){
    // Check the vector lengths
    if(x.size() != y.size()){
        std::cout << "Vectors must have the same length" << std::endl;
        return 1; // Error code
    }
    // Create the file (or erase an existing one)
    fileName = "results/" + fileName + ".txt";
    const char *fileNameChar = fileName.c_str();
    FILE *fp = fopen(fileNameChar, "w");
    if(fp != NULL){
        // Write the x values
        for(int j=0 ; j < x[0].size() ; j++){
            for(int i=0 ; i < x.size() ; i++)
                fprintf(fp, "%f ", x[i][j]);
            fprintf(fp, "\n");
        }
        // Write the y values
        for(int i=0 ; i < y.size() ; i++)
            fprintf(fp, "%f ", y[i]);
        fclose(fp);
        return 0;
    }
    else{
        printf("Cannot create a .txt file\n");
        return 2; // Error code
    }
}
