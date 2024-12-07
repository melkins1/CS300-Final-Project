#ifndef _hamming_h
#define _hamming_h

#include <eigen3/Eigen/Dense>
#include <fstream>
#include <string>
#include <utility>

class Hamming {
 public:
  //! Encodes a set of 4-bit messages into 7-bit messages with 3 parity bits
  /*! Each message is encoded through the following process: Transpose the message to be a 4x1 matrix, take the product of this Hamming code's generator matrix and the transposed message, turn every integer in the product into a bit by taking modulo 2, and transpose the result. The encoded version of the message has parity bits in positions 1, 2, and 4. The first parity bit gets the parity of the sum of bits 1, 2, and 4; the second parity bit gets the parity of the sum of bits 1, 3, and 4; and the third parity bit gets the parity of the sum of bits 2, 3, and 4.
   * @param data A matrix where each row is a 4-bit message
   * @return A matrix in which each row is the encoded version of the original matrix's row
   */
  Eigen::Matrix<bool,Eigen::Dynamic,7> encode(const Eigen::Matrix<bool,Eigen::Dynamic,4> &data) const;
  //! Decodes a set of 7-bit messages with 3 parity bits into the original 4-bit messages, with the ability to detect and correct a bitflip error if no more than one bit was flipped
  /*! Each message is decoded through the following process: Transpose the message to be a 7x1 matrix, take the product of the Hamming code's parity check matrix and the transposed message to get a syndrome vector, turn every integer in the syndrome vector into a bit by taking modulo 2, and get a wrong bit index by reversing the syndrome vector, converting it from binary to decimal, and subtracting 1. The parity check and generator matrices are set up such that this process will return -1 if no bit was flipped, and if a bit was flipped, it will return the index where the flip occurred instead. The program then corrects the bitflip and removes the parity bits, resulting in the original unencoded message. This process will fail if more than one bit was flipped.
   * @param data A matrix where each row is a 7-bit message with 3 parity bits
   * @return A matrix in which each row is the decoded version of the original matrix's row with any single bit-flip error corrected
   */
  Eigen::Matrix<bool,Eigen::Dynamic,4> correct(Eigen::Matrix<bool,Eigen::Dynamic,7> data) const;
  //! Gets a matrix of bits from a file
  /*! Each line of the file must be a set of either 4 or 7 bits, not separated by any spaces. If the file consists of 4-bit messages, you can encode them into 7-bit messages with 3 parity bits. If the file consists of 7-bit messages, you can decode them into 4-bit messages that use the original message's parity bits to correct any error that resulted from a single bit flip.
   * @param file An ifstream of the file to be read from
   * @param encoding A bool that is true if the matrix is to be encoded and false if the matrix is to be decoded
   * @return A pair consisting of an error message (empty string if no error occurred) and the Eigen matrix representation of the file's contents (incomplete if an error occurred)
   */
  std::pair<std::string,Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic>> getMatrixFromFile(std::ifstream &file, const bool &encoding) const;
 private:
  static const int pIndices[3][3];
  static const int pColumns[3];
  static const int nonPColumns[4];
  static const Eigen::Matrix<bool,7,4> generator;
  static const Eigen::Matrix<bool,3,7> parityCheck;
};

#endif
