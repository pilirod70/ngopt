/**
 * @file calcCoverMummer.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"

#include "OptionSet.h"
#include "Sequence.h"
#include "SequenceReader.h"
#include "SequenceWriter.h"
#include "Utils.h"
#include "Contig.h"
#include "BlatRecord.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <deque>

using namespace std;

char buf[MaxLine];
char line[MaxLine];
char name[MaxLine];
char ref_name[MaxLine];
vector<int64> recrod;
vector<int64> lengths;
vector<int64> valid_lengths;
set<string> valid_contigs;
map<string, int> dict;
int min_contig = 0;
double similar = 0.95;
double component_similar = 0.95;
vector<int> flags[10000];
int num_ref = 0;
map<string, int> species_id;
OptionSet option_set;
int num_species = 1;
vector<string> names;
vector<Sequence> contigs;

int main(int argc, char *argv[])
{
    option_set.AddOption("minContig", "", min_contig, "");
    option_set.AddOption("similar", "", similar, "");
    option_set.AddOption("component_similar", "", component_similar, "");
    option_set.AddOption("num", "n", num_species, "");
    option_set.ProcessOptions(argc, argv);

    if (argc < 2 + num_species)
    {
        cout << "usage validateComponent [--num n] component-file species-file1 species-file2" << endl;
        exit(1);
    }

    deque<int> component_size;

    FastaReader contig_reader(argv[1]);
    Sequence contig;
    string comment;
    while (contig_reader.Read(contig, comment))
    {
        names.push_back(comment);
        contigs.push_back(contig);
        int x, y;
        sscanf(comment.c_str() + 9, "%d_%d", &x, &y);
        
        if (x >= component_size.size())
            component_size.resize(x + 1);

        if ( contig.Size() >= min_contig)
            component_size[x] += contig.Size();
    }

    Sequence seq;
    for (int i = 0; i < num_species; ++i)
    {
        FastaReader reader(argv[2 + i]);
        
        while (reader.Read(seq, comment))
            species_id[comment] = i;
    }


    fgets(line, MaxLine, stdin);
    fgets(line, MaxLine, stdin);
    fgets(line, MaxLine, stdin);
    fgets(line, MaxLine, stdin);
    fgets(line, MaxLine, stdin);

    deque<deque<int> > votes(component_size.size());
    for (unsigned i = 0; i < votes.size(); ++i)
        votes[i].resize(num_species);

    int index = 0;
    map<string, map<int, int> > dict;
    while (fgets(line, MaxLine, stdin) != NULL)
    {
        BlatRecord record;
        record.Parse(line);

        int x, y;
        sscanf(record.query_name.c_str() + 9, "%d_%d", &x, &y);
        if (record.match_count > similar * record.query_length && record.query_length >= min_contig)
        {
            //cout << record.query_name << endl;
            //cout << x << " " << y << " " << species_id[record.ref_name] << " " << record.match_count << endl;
            if (dict[record.query_name][species_id[record.ref_name]] == 0)
            {
                if (record.match_count * 1.0 / record.query_length < 0.90)
                {
                    //cout << record.match_count << " " << record.query_length << " " << record.match_count * 100.0 / record.query_length << endl;
//                    for (unsigned i = 0; i < record.blocks.size(); ++i)
//                        cout << record.blocks[i].query_from << " " << record.blocks[i].size << endl;
                }
                votes[x][species_id[record.ref_name]] += record.match_count;
                dict[record.query_name][species_id[record.ref_name]] = 1;
            }
        }
    }

    int good = 0;
    int good_length = 0;
    int total_length = 0;
    for (unsigned i = 0; i < votes.size(); ++i)
    {
        if (*max_element(votes[i].begin(), votes[i].end()) > component_similar * component_size[i])
        {
            ++good;
            good_length += component_size[i];
            //good_length += *max_element(votes[i].begin(), votes[i].end());
        }

        //cout << votes[i][0] << " " << votes[i][1] << endl;

        total_length += component_size[i];
    }

    cout << FormatString("good %d/%d good_length %d/%d", good, component_size.size(), good_length, total_length) << endl;

    return 0;
}

