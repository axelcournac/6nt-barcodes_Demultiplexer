#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <vector>
#include <sys/stat.h>

using namespace std;

/*Barcode objects implement barcodes from Bioo NEXTflex 6nt barcodes (48 barcodes kit) used for reads demultiplexing.
 * Attributes: 
 *  m_bcName, the name of the bar code
 *  m_bcSeq, the seq of the barcode
 * Methods:
 *  getBCName(), getBCSeq: the getters*/
class Barcode {
    public:
        // Constructor
        Barcode(string bcName, string bcSeq):
            m_bcName(bcName), m_bcSeq(bcSeq) {}
        
        // Getters
        string getBCName() const {
            return m_bcName;
        }        
        string getBCSeq() const {
            return m_bcSeq;
        }
            
    private:
        string m_bcName;
        string m_bcSeq;
};


/*Barcodes object is mainly a vector of all the Bioo NEXTflex 6nt barcodes.
 * Attributes: 
 *  m_barcodes, a vector of barcodes
 * Methods:
 *  Barcodes(string barcodesFile): the constructor. Takes the barcode.txt file in entry.
 *                                 Must contain 1 barcode/line whith the comma-separated
 *                                 name and sequence of the barcode. The name must have 4 chars
 *                                 the sequence must have 6 characters.
 *  getBCs(): the getter of the vector*/
class Barcodes {
    public:
        // Constructor
        Barcodes(string barcodesFile): m_barcodesFile(barcodesFile) {
            initializeBarcodes();  
        }

        // Getters
        vector<Barcode> getBCs() {
            return m_barcodes;
        }
        
        // Other
        void initializeBarcodes() {
            ifstream flow(m_barcodesFile);
            if(!flow) {
                cerr << "ERROR: Could'n find the \"barcode.txt\" file." << endl;
            }
            
            else {
                string line, name, seq;
                // parse each line of the barcode file
                while(getline(flow, line)) {
                    //get the name and sequence of the barcodes
                    name = line.substr(0, 4);
                    seq = line.substr(5, 6);
        
                    //create the corresponding Barcode
                    Barcode bc(name, seq);
                    //add them in the corresponding vector
                    m_barcodes.push_back(bc);
                }
            }
        }

    private:
        string m_barcodesFile;
        vector<Barcode> m_barcodes;
};

/*BDemultiplexer object is the object which will compute the actual demultiplexing of the fastq file given as argument.
 * Attributes:
 *  m_inFile, (string) the name of the input file
 *  m_end, (string) indicates if the end1 file or the end2 file will be demultiplexed
 *  m_barcodes, (vector<Barcode>) a vector of Barcode object built from the barcode.txt file.
 *  m_counts, (vector<int>) the vector keeping track of the number of reads in every file.
 * Methods:
 *  Demultiplexer(string inFile, vector<Barcode>& barcodes, string end)the constructor.
 *  getBCs(): the getter of the input file name
 *  demultiplex(): the actual demultiplexer function. Open one ofstream per barcode + 1 for undetermined reads.
 *                 Gets each read and their barcode and writes them in the corresponding file thanks to the correct ofstream.
 *                 If the read's barcode is not in the vector, writes it in the file for undetermined reads.
 * */

class Demultiplexer {
    public:
        // Constructor
        Demultiplexer(string inFile, vector<Barcode>& barcodes, string end): 
            m_inFile(inFile), m_barcodes(barcodes), m_end(end) {
                vector<int> m_counts;
        }

        // Getter
        string getInfile() const {
            return m_inFile;
        }

        // Other
        //~ bool reportRead(Barcode barcode, ofstream ofStream, string end, string l1, string l2, string l3, string l4) {
            //~ ofStream << l1 << endl << l2 << endl << l3 << endl << l4 << endl;
            //~ return true;
        //~ }
        void report(string m_inFile, long long int ct, vector<Barcode> m_barcodes, vector<long long int> m_counts) {
            cout << "Reads of " << m_inFile << " demultiplexed:" << endl;
            cout << "Total: " << ct << endl;
            cout << "Name Seq Nb-of-reads" << endl;
            for(size_t i(0); i < m_barcodes.size(); i++) {
                cout << m_barcodes[i].getBCName() << " " << m_barcodes[i].getBCSeq() << " " << m_counts[i] << endl;
            }
            cout << "Undetermined XXXXXX " << m_counts[m_counts.size()-1] << endl;
        }

        long long int demultiplex() {
                //Ouverture d'un fichier en lecture
                ifstream inFlow(m_inFile);

                if(inFlow) {
                    // Declare variables
                    long long int ct = 0; // count the number of reads processed
                    //~ Barcode undeterminedBarcode("Undetermined", "......"); //a barcode used to export the reads whithout correct barcode in the Undetermined file
                    //ofstream ofStream;
                    string l1, l2, l3, l4;
                    bool reported;
                    int m_countsSize;

                    //Open 48 ofstreams to write into them
                    vector<FILE*> ofStreams;
                    for(auto bc: m_barcodes) {
                        string fileName = "Demultiplexed_Reads/" + bc.getBCName() + "." +  m_end + "fastq";
                        FILE* ofs;
                        ofs = fopen(fileName.c_str(), "w");
                        ofStreams.push_back(ofs);
                        m_counts.push_back(0);
                    }
                    string undeterminedOfStreamName = "Demultiplexed_Reads/Undetermined." +  m_end + ".fastq";
                    FILE* undeterminedOfStream;
                    undeterminedOfStream = fopen(undeterminedOfStreamName.c_str(), "w");
                    ofStreams.push_back(undeterminedOfStream);
                    m_counts.push_back(0);
                    m_countsSize = m_counts.size();

                    // Read the file 4 lines at a time and demultiplex as long as there are lines to read
                    while(getline(inFlow, l1) && getline(inFlow, l2) && getline(inFlow, l3) && getline(inFlow, l4)) {
                        
                        //get the current read's barcode value
                        string read_barcode = l1.substr(l1.find(" ")+7, 6);
                        reported = false;
                        
                        //check if it corresponds to a NEXTfera barcode
                        for(size_t i(0); i < m_barcodes.size(); i++) {
                            //if it does, write the read in the corresponding file and return true
                            if(read_barcode == m_barcodes[i].getBCSeq()) {
                                string toReport(l1 + "\n" + l2 + "\n" + l3 + "\n" + l4 + "\n");
                                int written = fwrite(toReport.c_str(), 1, toReport.length(), ofStreams[i]);
                                if (written != toReport.length()) {
                                    cerr << "Write error (probably end of disk)" << endl;
                                    exit(EXIT_FAILURE);
                                }
                                m_counts[i] = m_counts[i] + 1;
                                reported = true;
                                break;
                            }
                        }
                        //if not, write it in the "Undetermined" file
                        if(!reported) {
                            string toReport(l1 + "\n" + l2 + "\n" + l3 + "\n" + l4 + "\n");
                            int written = fwrite(toReport.c_str(), 1, toReport.length(), ofStreams[m_countsSize-1]);
                            if (written != toReport.length()) {
                                    cerr << "Write error (probably end of disk)" << endl;
                                    exit(EXIT_FAILURE);
                            }
                            m_counts[m_countsSize-1] = m_counts[m_countsSize-1] + 1;
                        }
                        
                        ct += 1;
                        if(ct%1000000 == 0) {
                            cout << ct/1000000 << " million reads\r" << flush;
                        }
                    }
                
                    report(m_inFile, ct, m_barcodes, m_counts);

                }
                
                else {
                    cerr << "ERROR: The file to demultiplex can't be opened." << endl;
                }
            
            return 0;
    }
        
    private:
        string m_inFile;
        string m_end;
        vector<Barcode> m_barcodes;
        vector<long long int> m_counts;
};    


int main(int argc, char* argv[]) {
    
    // Check that 2 files have been given as arguments
    if(argc != 3){
        cerr << "ERROR: wrong number of arguments. Please provide 1 or 2 fastq files to demultiplex." << endl;
        return 1;
    }
    
    
    //Create a new "Demultiplexed_Reads" folder
    system("rm -rf Demultiplexed_Reads");
    if(mkdir("Demultiplexed_Reads", 0755) != 0) {
        cerr << "ERROR: Couldn't create the Demultiplexed folder." << endl;
        return 1;
    }
            
    else {
    
        // Get the files to demultiplex
        string f1 = argv[1];
        string f2 = argv[2];
        
        // Generate the vectors of barcodes
        Barcodes bcs1("barcodes.txt", "end1");
        vector<Barcode> barcodes1(bcs1.getBCs());
        Barcodes bcs2("barcodes.txt", "end2");
        vector<Barcode> barcodes2(bcs2.getBCs());
        
        // Demultiplex the reads of the end1
        Demultiplexer dm1(f1, barcodes1, "end1");
        dm1.demultiplex();
        //cout << "\nReads of end1 demultiplexed" << endl;
        
        // Demultiplex the reads of the end2
        if(f2 != "-") {
            Demultiplexer dm2(f2, barcodes2, "end2");
            dm2.demultiplex();
            //cout << "\nReads of end2 demultiplexed" << endl;
        }
        cout << endl << "Files saved in Demultiplexed_Reads directory." << endl;

        return 0;
    }
}

