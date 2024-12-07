#include "Hamming.h"
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>

Eigen::Matrix<bool,Eigen::Dynamic,7> Hamming::encode(const Eigen::Matrix<bool,Eigen::Dynamic,4> &data) const {
  const int numRows = data.rows();
  Eigen::Matrix<bool,Eigen::Dynamic,7> encoded;
  encoded.resize(numRows,7);
  for (int r=0; r<numRows; r++) {
    Eigen::Matrix<bool,4,1> transposedRow = data.row(r).transpose();
    Eigen::Matrix<int,7,1> newRowNotBits = generator.cast<int>()*transposedRow.cast<int>();
    Eigen::Matrix<bool,7,1> newRowUntransposed;
    for (int i=0; i<7; i++) {
      newRowUntransposed(i,0) = newRowNotBits(i,0)%2;
    }
    encoded.row(r) = newRowUntransposed.transpose();
  }
  return encoded;
}

Eigen::Matrix<bool,Eigen::Dynamic,4> Hamming::correct(Eigen::Matrix<bool,Eigen::Dynamic,7> data) const {
  const int numRows = data.rows();
  Eigen::Matrix<bool,Eigen::Dynamic,4> corrected;
  corrected.resize(numRows,4);
  // correct any flips that have occurred
  for (int r=0; r<numRows; r++) {
    Eigen::Matrix<bool,7,1> transposedRow = data.row(r).transpose();
    Eigen::Matrix<int,3,1> syndromeNotBits = parityCheck.cast<int>()*transposedRow.cast<int>();
    Eigen::Matrix<bool,7,1> syndrome;
    for (int i=0; i<3; i++) {
      syndrome(i,0) = syndromeNotBits(i,0)%2;
    }
    int wrongBitIndex = syndrome(0,0)+syndrome(1,0)*2+syndrome(2,0)*4-1;
    if (wrongBitIndex!=-1) {
      data(r,wrongBitIndex) = !data(r,wrongBitIndex);
    }
    for (int i=0; i<4; i++) {
      corrected(r,i) = data(r,nonPColumns[i]);
    }
  }
  return corrected; 
}

std::pair<std::string,Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>> Hamming::getMatrixFromFile(std::ifstream &file, const bool &encoding) const {
  std::pair<std::string,Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>> returnPair;
  std::string errorMessage = "";
  Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> messages;
  if (!file.good()) {
    errorMessage = "Error: file does not appear to exist";
    file.close();
    return std::make_pair(errorMessage,messages);
  }
  int numCols = encoding ? 4 : 7;
  int currentRow = 0;
  std::string line;
  while (std::getline(file,line)) {
    line.erase(std::remove(line.begin(),line.end(),'\r'),line.end());
    if (line.length() != numCols) {
      errorMessage = "Error: line "+std::to_string(currentRow+1)+" is of length "+std::to_string(line.length())+", not "+std::to_string(numCols);
      file.close();
      return std::make_pair(errorMessage,messages);
    }
    messages.conservativeResize(currentRow+1,numCols);
    for (int i=0; i<line.length(); i++) {
      char c = line[i];
      if (c!='0' && c!='1') {
	errorMessage = "Error: invalid character "+std::string(1,c)+" at line "+std::to_string(currentRow+1);
	file.close();
	return std::make_pair(errorMessage,messages);
      }
      else {
	messages(currentRow,i) = c=='1';
      }
    }
    currentRow += 1;
  }
  file.close();
  return std::make_pair(errorMessage,messages);
} 

const int Hamming::pIndices[3][3] = {{0,1,3},{0,2,3},{1,2,3}};
const int Hamming::pColumns[3] = {0,1,3};
const int Hamming::nonPColumns[4] = {2,4,5,6};
const Eigen::Matrix<bool,7,4> Hamming::generator((Eigen::Matrix<bool,7,4>() <<
  true, true, false, true,
  true, false, true, true,
  true, false, false, false,
  false, true, true, true,
  false, true, false, false,
  false, false, true, false,
  false, false, false, true).finished());
const Eigen::Matrix<bool,3,7> Hamming::parityCheck((Eigen::Matrix<bool,3,7>() <<
  true, false, true, false, true, false, true,
  false, true, true, false, false, true, true,
  false, false, false, true, true, true, true).finished());

std::string matrixRowToString(Eigen::Matrix<bool,1,Eigen::Dynamic> row) {
  std::string strToReturn = "";
  for (int i=0; i<row.cols(); i++) {
    char c = row(0,i) ? '1' : '0';
    strToReturn += c;
  }
  return strToReturn;
}

void unitTests() {
  // all tests should show true
  std::pair<std::string,Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>> testPair;
  Hamming hamming;
  // test rejection of files that don't exist
  std::ifstream testFile1("testFiles/thisFileDoesNotExist.txt");
  testPair = hamming.getMatrixFromFile(testFile1,true);
  std::cout << (testPair.first.find("exist") != std::string::npos) << std::endl;
  // test rejection of files with invalid characters
  std::ifstream testFile2("testFiles/invalidCharTest.txt");
  testPair = hamming.getMatrixFromFile(testFile2,true);
  std::cout << (testPair.first.find("character") != std::string::npos) << std::endl;
  // test rejection of files with improper length for procedure
  std::ifstream testFile3("testFiles/wrongLengthToEncode.txt");
  testPair = hamming.getMatrixFromFile(testFile3,true);
  std::cout << (testPair.first.find("length") != std::string::npos) << std::endl;
  std::ifstream testFile4("testFiles/wrongLengthToDecode.txt");
  testPair = hamming.getMatrixFromFile(testFile4,false);
  std::cout << (testPair.first.find("length") != std::string::npos) << std::endl;
  // test acceptance of files that should not have issues
  std::ifstream testFile5("testFiles/shouldEncodeProperly.txt");
  testPair = hamming.getMatrixFromFile(testFile5,true);
  std::cout << (testPair.first.length()==0) << std::endl;
  std::ifstream testFile6("testFiles/shouldDecodeProperly.txt");
  testPair = hamming.getMatrixFromFile(testFile6,false);
  std::cout << (testPair.first.length()==0) << std::endl;
  // test various message encoding
  std::ifstream testFile7Through9("testFiles/messagesToEncode.txt");
  testPair = hamming.getMatrixFromFile(testFile7Through9,true);
  Eigen::Matrix<bool,Eigen::Dynamic,7> testEncodedMatrix = hamming.encode(testPair.second);
  std::cout << (matrixRowToString(testEncodedMatrix.row(0))=="1011010") << std::endl;
  std::cout << (matrixRowToString(testEncodedMatrix.row(1))=="0010110") << std::endl;
  std::cout << (matrixRowToString(testEncodedMatrix.row(2))=="1101001") << std::endl;
  // test various message decoding (no errors present)
  std::ifstream testFile10Through12("testFiles/messagesToDecodeNoErrors.txt");
  testPair = hamming.getMatrixFromFile(testFile10Through12,false);
  Eigen::Matrix<bool,Eigen::Dynamic,4> testDecodedMatrix = hamming.correct(testPair.second);
  std::cout << (matrixRowToString(testDecodedMatrix.row(0))=="1010") << std::endl;
  std::cout << (matrixRowToString(testDecodedMatrix.row(1))=="1110") << std::endl;
  std::cout << (matrixRowToString(testDecodedMatrix.row(2))=="0001") << std::endl;
  std::ifstream testFile13Through19("testFiles/messagesToDecodeWithErrors.txt");
  testPair = hamming.getMatrixFromFile(testFile13Through19,false);
  Eigen::Matrix<bool,Eigen::Dynamic,4> testDecodedMatrix2 = hamming.correct(testPair.second);
  for (int i=0; i<7; i++) {
    std::cout << (matrixRowToString(testDecodedMatrix2.row(i))=="1010") << std::endl;
  }
}

int main() {
  //unitTests();
  std::cout << "Create a text file with lines of binary messages separated by line breaks. The individual bits should not have any spaces between them. If you are encoding messages, there should be 4 bits per message; if you are decoding messages, there should be 7 bits per message." << std::endl;
  std::cout << "Enter the path to this file (e.g. testFiles/exampleFile.txt)" << std::endl;
  std::string filename;
  std::getline(std::cin,filename);
  std::cout << "Enter e if you are encoding. Enter d if you are decoding" << std::endl;
  std::string modeStr;
  std::getline(std::cin,modeStr);
  while (modeStr[0]!='d'&&modeStr[0]!='e') {
    std::cout << "Invalid mode; please type e or d" << std::endl;
    std::getline(std::cin,modeStr);
  }
  bool mode = modeStr[0]=='e';
  Hamming hamming;
  std::ifstream file(filename);
  std::pair<std::string,Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>> outputPair = hamming.getMatrixFromFile(file,mode);
  if (outputPair.first.length()!=0) {
    std::cout << outputPair.first << std::endl;
  }
  else {
    if (mode) {
      std::cout << hamming.encode(outputPair.second) << std::endl;
    }
    else {
      std::cout << hamming.correct(outputPair.second) << std::endl;
    }
  }
  std::cout << "Restart the program to encode or decode another file" << std::endl;
  return 0;
}
