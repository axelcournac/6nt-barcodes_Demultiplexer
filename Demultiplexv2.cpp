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
 *                                 name and sequence of the barcode. The name must have 8 chars
 *                                 the sequence must have 17 characters.
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
                    name = line.substr(0, 8);
                    seq = line.substr(9, 17);
        
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
        Demultiplexer(string inFile1, string inFile2, vector<Barcode>& barcodes):
            m_inFile1(inFile1), m_inFile2(inFile2), m_barcodes(barcodes) {
                vector<int> m_counts1, m_counts2;
        }

        // Other 
        /*demultiplex(): the backbone of demultiplexing pipeline*/      
        int demultiplex() {
            //if f2 == "None", demultiplex file f1 as a single_end
            if(m_inFile2 == "None") {
                ifstream inFlow(m_inFile1); //open m_inFile1
                if(inFlow) {
                    // Declare variables
                    int m_countsSize; // a vector keeping track of how many reads are associated to each barcode
                    string end("");
                    //Open 48 ofstreams to write into them and put them in openOfStreams
                    vector<FILE*> ofStreams = openOfStreams(end);
                    m_countsSize = m_counts1.size();
                    // Read the file 4 lines at a time and demultiplex as long as there are lines to read
                    readAndSort(inFlow, ofStreams, m_countsSize);
                }
                else {
                    cerr << "ERROR: The file to demultiplex can't be opened." << endl;
                }
            }
            
            // else, demultiplex files f1 and f2 as pair-end files
            else {
                ifstream inFlow1(m_inFile1); //open m_inFile1
                ifstream inFlow2(m_inFile2); //open m_inFile2
                if(inFlow1 && inFlow2) {
                    // Declare variables
                    int m_countsSize;
                    string end(".end1");
                    //Open 48 ofstreams to write into them and put them in openOfStreams
                    vector<FILE*> ofStreams1 = openOfStreams(end);
                    end = ".end2";
                    vector<FILE*> ofStreams2 = openOfStreams(end);
                    m_countsSize = m_counts1.size();
                    // Read the file 4 lines at a time and demultiplex as long as there are lines to read
                    readAndSort(inFlow1, inFlow2, ofStreams1, ofStreams2, m_countsSize);
                }
            }
            return 0;
        }
        
        /*openOfStreams: create and returns a vector of streams to files corresponding to all the barcodes*/
        vector<FILE*> openOfStreams(string end) {                   
            //open the streams corresponding to each barcode
            vector<FILE*> ofStreams;
            for(auto bc: m_barcodes) {
                string fileName = "Demultiplexed_Reads/" + bc.getBCName() + end + ".fastq";
                FILE* ofs;
                ofs = fopen(fileName.c_str(), "w");
                ofStreams.push_back(ofs);
                if(end == "" || end == ".end1") {
                    m_counts1.push_back(0);
                }
                else{
                    m_counts2.push_back(0);
                }
            }
            //open a special stream for undertermined reads
            string undeterminedOfStreamName = "Demultiplexed_Reads/Undetermined." +  end + ".fastq";
            FILE* undeterminedOfStream;
            undeterminedOfStream = fopen(undeterminedOfStreamName.c_str(), "w");
            ofStreams.push_back(undeterminedOfStream);
            if(end == "" || end == ".end1") {
                m_counts1.push_back(0);
            }
            else{
                m_counts2.push_back(0);
            }

            return ofStreams;
        }
        
        /**/
        void readAndSort(ifstream& inFlow, vector<FILE*>& ofStreams, int m_countsSize) { 
            int written;
            long long int ct(0);
            string l1, l2, l3, l4, toReport, read_barcode;
            bool reported;
            while(getline(inFlow, l1) && getline(inFlow, l2) && getline(inFlow, l3) && getline(inFlow, l4)) {
                //get the current read's barcode value
                read_barcode = l1.substr(l1.find(" ")+7, 17);
                reported = false;

                //check if it corresponds to a NEXTfera barcode
                for(size_t i(0); i < m_barcodes.size(); i++) {
                    //if it does, write the read in the corresponding file and return true
                    if(read_barcode == m_barcodes[i].getBCSeq()) {
                        toReport = l1 + "\n" + l2 + "\n" + l3 + "\n" + l4 + "\n";
                        written = fwrite(toReport.c_str(), 1, toReport.length(), ofStreams[i]);
                        if (written != toReport.length()) {
                            cerr << "Write error (probably end of disk)" << endl;
                            exit(EXIT_FAILURE);
                        }
                        m_counts1[i] = m_counts1[i] + 1;
                        reported = true;
                        break;
                    }
                }
                
                //if not, write it in the "Undetermined" file
                if(!reported) {
                    toReport = l1 + "\n" + l2 + "\n" + l3 + "\n" + l4 + "\n";
                    int written = fwrite(toReport.c_str(), 1, toReport.length(), ofStreams[m_countsSize-1]);
                    if (written != toReport.length()) {
                        cerr << "Write error (probably end of disk)" << endl;
                        exit(EXIT_FAILURE);
                    }
                    m_counts1[m_countsSize-1] = m_counts1[m_countsSize-1] + 1;
                }
                            
                ct += 1;
                if(ct%1000000 == 0) {
                    cout << ct/1000000 << " million reads\r" << flush;
                }
            }
                    
            report(m_inFile1, ct, m_barcodes, m_counts1);
        }
        
        void readAndSort(ifstream& inFlow1, ifstream& inFlow2, vector<FILE*>& ofStreams1, vector<FILE*>& ofStreams2, int m_countsSize) { 
            int written1, written2;
            long long int ct(0);
            string l11, l12, l13, l14, l21, l22, l23, l24, toReport1, toReport2, read_barcode1, read_barcode2;
            bool reported, notSameBarcode(false);
            while(getline(inFlow1, l11) && getline(inFlow1, l12) && getline(inFlow1, l13) && getline(inFlow1, l14) && \
                getline(inFlow2, l21) && getline(inFlow2, l22) && getline(inFlow2, l23) && getline(inFlow2, l24)) {
                //get the current read's barcode value
                read_barcode1 = l11.substr(l11.find(" ")+7, 17);
                read_barcode2 = l21.substr(l21.find(" ")+7, 17);
                reported = false;
                //if barcode1 and barcode2 are the same, check if they correspond to a known barcode
                if(read_barcode1 == read_barcode2) {
                    for(size_t i(0); i < m_barcodes.size(); i++) {
                        //if it does, write the read in the corresponding file and return true
                        if(read_barcode1 == m_barcodes[i].getBCSeq()) {
                            toReport1 = l11 + "\n" + l12 + "\n" + l13 + "\n" + l14 + "\n";
                            toReport2 = l21 + "\n" + l22 + "\n" + l23 + "\n" + l24 + "\n";
                            written1 = fwrite(toReport1.c_str(), 1, toReport1.length(), ofStreams1[i]);
                            written2 = fwrite(toReport2.c_str(), 1, toReport2.length(), ofStreams2[i]);
                            if(written1 != toReport1.length() || written2 != toReport2.length()) {
                                cerr << "Write error (probably end of disk)" << endl;
                                exit(EXIT_FAILURE);
                            }
                            m_counts1[i] = m_counts1[i] + 1;
                            m_counts2[i] = m_counts2[i] + 1;
                            reported = true;
                            break;
                            
                        }
                    }
                }
                
                //if not, keep it in memory
                else {
                    notSameBarcode = true;
                    reported = false;
                }
                
                //if the reads have not been reported for any reason, write them in the "Undetermined" file
                if(!reported){
                    toReport1 = l11 + "\n" + l12 + "\n" + l13 + "\n" + l14 + "\n";
                    toReport2 = l21 + "\n" + l22 + "\n" + l23 + "\n" + l24 + "\n";
                    written1 = fwrite(toReport1.c_str(), 1, toReport1.length(), ofStreams1[m_countsSize-1]);
                    written2 = fwrite(toReport2.c_str(), 1, toReport2.length(), ofStreams2[m_countsSize-1]);
                    if(written1 != toReport1.length() || written2 != toReport2.length()) {
                        cerr << "Write error (probably end of disk)" << endl;
                        exit(EXIT_FAILURE);
                    }
                    m_counts1[m_countsSize-1] = m_counts1[m_countsSize-1] + 1;
                    m_counts2[m_countsSize-1] = m_counts2[m_countsSize-1] + 1;
                }
                            
                ct += 1;
                if(ct%1000000 == 0) {
                    cout << ct/1000000 << " million reads\r" << flush;
                }
            }
            
            if(notSameBarcode) {
                cout << "Some pair-end reads didn't have the same barcode !" << endl;
            }
            report(m_inFile1, ct, m_barcodes, m_counts1);
            report(m_inFile2, ct, m_barcodes, m_counts2);
        }
        
        /*display a short report of the demultipeling*/
        void report(string m_inFile, long long int ct, vector<Barcode> m_barcodes, vector<long long int> m_counts) {
            cout << "Reads of " << m_inFile << " demultiplexed:" << endl;
            cout << "Total: " << ct << endl;
            cout << "Name Seq Nb-of-reads" << endl;
            for(size_t i(0); i < m_barcodes.size(); i++) {
                cout << m_barcodes[i].getBCName() << " " << m_barcodes[i].getBCSeq() << " " << m_counts[i] << endl;
            }
            cout << "Undetermined XXXXXX " << m_counts[m_counts.size()-1] << endl;
        }


    private:
        string m_inFile1;
        string m_inFile2;
        vector<Barcode> m_barcodes;
        vector<long long int> m_counts1;
        vector<long long int> m_counts2;
};    






int main(int argc, char* argv[]) {
    
    // Check that 2 files have been given as arguments
    if(argc > 3 || argc <= 1 ){
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
        string f1(argv[1]);
        string f2;
        if(argc == 3) {
            f2 = argv[2];
        }
        else {
            f2 = "None";
        }
        
        // Generate the vector of barcodes
        Barcodes bcs("names_barcodes_2019.txt");
        vector<Barcode> barcodes(bcs.getBCs());
        
        // Demultiplex the reads of the end1
        Demultiplexer dm(f1, f2, barcodes);
        dm.demultiplex();
        
        return 0;
    }
}

