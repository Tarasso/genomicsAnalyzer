#include<iostream>
#include<fstream>
#include<string>
#include<iomanip>
#include<sstream>
#include<vector>
#include<cmath>
#include<algorithm>

using namespace std;

int size;

// Struct thatr stores the data for each gene
struct GENE
{
  string description;
  int start;
  int end;
  char strand;
  int numProteins;
  string unNeeded1;
  string unNeeded2;
  string name;
  string orfid;
};

// Method to split lines on a delimiter
vector<string> split(string s, char c)
{
  istringstream ss(s);
  string token;
  vector<string> tokens;
  while(getline(ss,token,c))
    tokens.push_back(token);
  return tokens;
}

// Reads Fasta file and returns a string of all base pairs and prints header of file
string readFasta(char* fileName)
{
  string genome;
  ifstream file;
  file.open(fileName);
  if(!file.is_open()) cout << "Could not open fasta file" << endl;
  string info;
  getline(file,info);
  cout << "Reading: " << info.substr(1) << endl;
  string line;
  while(getline(file,line))
  {
    genome += line;
  }
  return genome;
}

// Calculates the percent of G's and C's in a genome
float GCcontent(string genome)
{
  int count;
  for(int i = 0; i < genome.length(); i++)
  {
    if(genome[i] == 'G' || genome[i] == 'C')
      count++;
  }
  return 100*((float)count/genome.length());
}

// Given a position and vector of all the genes, outputs the closest gene before the position
string findPrevGene(int position, vector<GENE> genes)
{
  GENE targetGene;
  for(int i = 1; i < genes.size(); i++)
  {
    if(genes.at(i).end>position)
    {
      targetGene = genes.at(i-1);
      break;
    }
  }
  string answer;
  answer = " " + targetGene.orfid + "/" + targetGene.name;
  return answer;
}

// Find the longest consecutive set of base pairs
string homopolymer(string genome, vector<GENE> genes, bool geneIndex = false)
{
  int largestString = -1;
  char base;
  int position;
  for(int i = 0; i < genome.length(); i++)
  {
    bool nextSame = true;
    int count = 0;
    while(nextSame)
    {
      count++;
      if(genome[i]==genome[i+count])
        nextSame = true;
      else
        nextSame = false;
    }
    if(count>largestString)
    {
      largestString = count;
      base = genome[i];
      position = i;
    }
  }
  string answer;
  for(int i = 0; i < largestString; i++)
  {
    answer += base;
  }
  answer = answer + " (len=" + to_string(largestString) + " bp) at coord " + to_string(position);
  if(geneIndex) answer = answer + ", after" + findPrevGene(position, genes); // Appends a gene reference location is geneIndex flag is true
  return answer;
}

// Given a Prot file, returns a vector of GENE structs including information about each gene
vector<GENE> readProt(char* filename)
{
  ifstream file;
  file.open(filename);
  if(!file.is_open()) cout << "Could not open prot_table file" << endl;
  string line;
  vector<GENE> genes;
  while(getline(file,line))
  {
    vector<string> vars = split(line,'\t');
    GENE tempGene;
    tempGene.description = vars[0];
    tempGene.start = atoi(vars[1].c_str());
    tempGene.end = atoi(vars[2].c_str());
    tempGene.strand = vars[3][0];
    tempGene.numProteins = atoi(vars[4].c_str());
    tempGene.unNeeded1 = vars[5];
    tempGene.unNeeded2 = vars[6];
    tempGene.name = vars[7];
    tempGene.orfid = vars[8];
    genes.push_back(tempGene);
  }
  return genes;
}

// Returns an integer vector including the smallest and largest gene in a vector of GENE
vector<int> smallLargeGene(vector<GENE> genes)
{
  int smallest = 999999;
  int largest = -999999;
  for(int i = 0; i < genes.size(); i++)
  {
    int size = genes.at(i).end - genes.at(i).start + 1;
    if(size < smallest)
      smallest = size;
    if(size > largest)
      largest = size;
  }
  vector<int> nums;
  nums.push_back(smallest);
  nums.push_back(largest);
  return nums;
}

// Calculates the mean and standard deviation of all the genes returning them in a vector of floats
vector<float> meanAndStd(vector<GENE> genes)
{
  int count = 0;
  float sum = 0;
  float squaredSum = 0;

  for(int i = 0; i < genes.size(); i++)
  {
    int size = genes.at(i).end - genes.at(i).start + 1;
    count++;
    sum += size;
    squaredSum += size * size;
  }

  float mean = sum/count;
  float stdev = sqrt((squaredSum/count)-mean*mean);

  vector<float> answer;
  answer.push_back(mean);
  answer.push_back(stdev);

  return answer;

}

// given a vector of GENE's, returns a vector of floats containing the percent coded, min, max, and mean of intergenic regions
vector<float> intergenicStats(vector<GENE> genes, string& minLocation, string& maxLocation, string& minName, string& maxName)
{
  int totalIntergenicIncludingNeg = 0;
  int maxIntergenic = -999999;
  int minIntergenic = 999999;
  for(int i = 0; i < genes.size()-1; i++)
  {
    int difference = genes.at(i+1).start - genes.at(i).end;
    totalIntergenicIncludingNeg += difference;
    if(difference > maxIntergenic)
    {
      maxIntergenic = difference + 1;
      maxLocation = genes.at(i).orfid;
      maxName = genes.at(i).name;
    }
    if(difference < minIntergenic)
    {
      minIntergenic = difference + 1;
      minLocation = genes.at(i).orfid;
      minName = genes.at(i).name;
    }
  }

  // also need to check distance between last gene and end of base pair sequence
  int lastDifference = size - genes.at(genes.size()-1).end;
  totalIntergenicIncludingNeg += lastDifference;
  if(lastDifference > maxIntergenic)
  {
    maxIntergenic = lastDifference;
    maxLocation = genes.at(genes.size()-1).orfid;
    maxName = genes.at(genes.size()-1).name;
  }
  if(lastDifference < minIntergenic)
  {
    minIntergenic = lastDifference;
    minLocation = genes.at(genes.size()-1).orfid;
    minName = genes.at(genes.size()-1).name;
  }

  float percentCoded = ((float)(size - totalIntergenicIncludingNeg) / size) * 100;

  float meanIntergenic = ((float) totalIntergenicIncludingNeg / genes.size());

  vector<float> answer;
  answer.push_back(percentCoded);
  answer.push_back(minIntergenic);
  answer.push_back(maxIntergenic);
  answer.push_back(meanIntergenic);

  return answer;

}

// helper methods that returns a genes index in the GENE vector given a gene ID or name
int findTargetGene(vector<GENE> genes, string nameOrID)
{
  for(int i = 0; i < genes.size(); i++)
  {
    if(genes.at(i).orfid == nameOrID || genes.at(i).name == nameOrID)
      return i;
  }
  return -1;
}


// method that prints a string formatted to so many character per line
void print_seq(string seq)
{
  int count = 0;
  int printCount = 1;
  for(int i = 0; i < seq.size(); i++)
  {
    if(i==0)
      cout << setw(4) << "1" << " ";
    cout << seq[i];
    count++;
    printCount++;
    if(count == 70) // 70 characters per line
    {
      cout << endl << setw(4) << printCount << " ";
      count = 0;
    }
  }
}

// returns a string containing a gene's base pairs given the genes start and end location and the fasta string
string getGeneSeq(int start, int end, string fastaGenome)
{
  string seq;
  for(int i = start; i <= end; i++)
  {
    seq += fastaGenome[i];
  }
  return seq;
}

// returns the reverse compliment of a given string of base pairs
string reverseCompliment(string seq)
{
  reverse(seq.begin(),seq.end());
  string newString;
  for(int i = 0; i < seq.size(); i++)
  {
    char base = seq[i];
    if (base=='A')
    {
      newString += 'T';
    }
    else if (base=='T')
    {
      newString += 'A';
    }
    else if (base=='G')
    {
      newString += 'C';
    }
    else if (base=='C')
    {
      newString += 'G';
    }
  }
  return newString;
}


// given a set of three base pairs as a string, returns the coresponding aminoacid
char codon (string s)
{
  if (s=="TTA") return 'L';
  if (s=="TTG") return 'L';
  if (s=="CTT") return 'L';
  if (s=="CTC") return 'L';
  if (s=="CTA") return 'L';
  if (s=="CTG") return 'L';
  if (s=="TGG") return 'W';
  if (s=="TAA") return '*';
  if (s=="TAG") return '*';
  if (s=="TGA") return '*';
  if (s=="ATG") return 'M';
  if (s=="TTT") return 'F';
  if (s=="TTC") return 'F';
  if (s=="TAT") return 'Y';
  if (s=="TAC") return 'Y';
  if (s=="TCT") return 'S';
  if (s=="TCC") return 'S';
  if (s=="TCA") return 'S';
  if (s=="TCG") return 'S';
  if (s=="AGT") return 'S';
  if (s=="AGC") return 'S';
  if (s=="CCT") return 'P';
  if (s=="CCC") return 'P';
  if (s=="CCA") return 'P';
  if (s=="CCG") return 'P';
  if (s=="TGT") return 'C';
  if (s=="TGC") return 'C';
  if (s=="CAT") return 'H';
  if (s=="CAC") return 'H';
  if (s=="CAA") return 'Q';
  if (s=="CAG") return 'Q';
  if (s=="AAT") return 'N';
  if (s=="AAC") return 'N';
  if (s=="CGT") return 'R';
  if (s=="CGC") return 'R';
  if (s=="CGA") return 'R';
  if (s=="CGG") return 'R';
  if (s=="AGA") return 'R';
  if (s=="AGG") return 'R';
  if (s=="ATT") return 'I';
  if (s=="ATC") return 'I';
  if (s=="ATA") return 'I';
  if (s=="AAA") return 'K';
  if (s=="AAG") return 'K';
  if (s=="GAT") return 'D';
  if (s=="GAC") return 'D';
  if (s=="GAA") return 'E';
  if (s=="GAG") return 'E';
  if (s=="ACT") return 'T';
  if (s=="ACC") return 'T';
  if (s=="ACA") return 'T';
  if (s=="ACG") return 'T';
  if (s=="GTT") return 'V';
  if (s=="GTC") return 'V';
  if (s=="GTA") return 'V';
  if (s=="GTG") return 'V';
  if (s=="GCT") return 'A';
  if (s=="GCC") return 'A';
  if (s=="GCA") return 'A';
  if (s=="GCG") return 'A';
  if (s=="GGT") return 'G';
  if (s=="GGC") return 'G';
  if (s=="GGA") return 'G';
  if (s=="GGG") return 'G';
  return '?';
}

// converts gene sequence to amino acids
string aminoAcids(string geneSeq)
{
  string ret;
  for(int i = 0; i < geneSeq.size(); i+=3)
  {
    string tempCodon;
    char c1 = geneSeq[i];
    char c2 = geneSeq[i+1];
    char c3 = geneSeq[i+2];
    tempCodon.append(1,c1);
    tempCodon.append(1,c2);
    tempCodon.append(1,c3);
    ret += codon(tempCodon);
  }
  ret[0] = 'M';
  return ret;
}

int main(int argc, char** argv)
{
  // reads files and sets the global variable size to the base pair length
  string fastaGenome = readFasta(argv[1]);
  vector<GENE> genes = readProt(argv[2]);
  size = fastaGenome.length();

  // all the method calls and print statements for the gene information
  cout << "Length = " << size << " bp" << endl;
  float GCpercent = GCcontent(fastaGenome);
  cout << "GC content = " << setprecision(3) << GCpercent << "%" << endl;
  string longestHomopolymer = homopolymer(fastaGenome, genes);
  cout << "longest homopolymer: " << longestHomopolymer << endl;
  cout << "num genes: " << genes.size() << endl;
  vector<int> minAndMax = smallLargeGene(genes);
  cout << "gene sizes: [" << minAndMax.at(0) << "," << minAndMax.at(1) << "]";
  vector<float> meanStd = meanAndStd(genes);
  cout << ", mean = " << fixed << setprecision(1) << meanStd.at(0) << " bp stdev = " << meanStd.at(1) << endl;

  // all the method calls and print statements for the intergenic information
  string minLocation;
  string maxLocation;
  string minName;
  string maxName;
  vector<float> intergenicData = intergenicStats(genes, minLocation, maxLocation, minName, maxName);
  cout << "coding fraction: " << fixed << setprecision(1) << intergenicData.at(0) << "%" << endl;
  cout << "mean size of intergenic regions: " << fixed << setprecision(1) << intergenicData.at(3) << " bp" << endl;
  cout << "largest intergenic: " << setprecision(0) << intergenicData.at(2) << " bp (after " << maxLocation << "/" << maxName << ")" << endl;
  cout << "smallest intergenic: " << intergenicData.at(1) << " bp (after " << minLocation << "/" << minName << ")" << endl;
  longestHomopolymer = homopolymer(fastaGenome,genes,true);
  cout << "longest homopolymer: " << longestHomopolymer << endl;
  cout << endl;

  // all the calls and print statements for the information regarding a specific gene
  if(argv[3] != NULL)
  {
    int targetGeneIndex = findTargetGene(genes, argv[3]);
    if(targetGeneIndex == -1)
      cout << "Gene " << argv[3] << " cannot be found." << endl;
    else
    {
      bool onNegStrand = false;
      int start = genes.at(targetGeneIndex).start - 1;
      int end = genes.at(targetGeneIndex).end - 1;
      string geneSeq = getGeneSeq(start,end,fastaGenome);

      if(genes.at(targetGeneIndex).strand == '-')
        onNegStrand = true;
      if(onNegStrand)
        geneSeq = reverseCompliment(geneSeq);

      print_seq(geneSeq);
      cout << endl << endl;
      string aminos = aminoAcids(geneSeq);
      print_seq(aminos);
    }
  }
  // error handling if gene cannot be found
  else
  {
    cout << "Could not print information about gene. No gene given to analyze." << endl;
  }
}
