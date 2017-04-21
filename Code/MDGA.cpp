#define _CRT_SECURE_NO_WARNINGS

#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdio>
#include<cmath>
#include<cstdlib>
#include<cctype>
#include<cstring>
#include<climits>
#include<iomanip>
#include<string>
#include<vector>
#include<stack>
#include<queue>
#include<deque>
#include<set>
#include<map>
#include<list>
#include<algorithm>
#include<utility>
#include<functional>

using namespace std;

#define psb(x) push_back(x)
#define psf(x) push_front(x)
#define ppb pop_back()
#define ppf pop_front()
#define pop pop()
#define front front()
#define back back()
#define bgn begin()
#define end end()
#define emp empty()
#define clr clear()
#define sz size()
#define sp setprecision
#define fx fixed
#define fst first
#define snd second

#define max3(x, y, z) max(x, max(y, z))
#define min3(x, y, z) min(x, min(y, z))
#define max4(w, x, y, z) max(w, max(x, max(y,z)))
#define min4(w, x, y, z) min(w, min(x, min(y, z)))

#define loop(x,r,n) for(x = r ; x <= n ; x++)
#define loopr(x, r, n) for (x = r; x >= n; x--)
#define test(t) for(int o = 1 ; o <= t ; o++)
#define printcs cout << "Case " << o << ": ";
#define nl cout << "\n"

#define INFMX 2139062143
#define INFMN -2139062144
#define pi acos(-1.0)
#define N 1000005

char base_pair[] = { 'A', 'T', 'G', 'C' };

struct WEIGHT_MATRIX
{
	vector<double> wm[4];
};

struct PATTERN
{
	string pat;
	double tfs;
	int complete_match;
	double percentage;
};

struct GENERATION_INFO
{
	int number_of_pattern,max_match;
	double max_tfs, min_tfs;
	vector<PATTERN> patterns;
	map<string, int> unique_pattern;
	map<string, vector<PATTERN> > best_matches;
	map<string, WEIGHT_MATRIX> weight_matrix_all;
	map<string, pair<string, string> > parents;
	map<string, pair<string, string> > children;
};

void print_patterns(vector<PATTERN> test);
double if_not_match(char base1, char base2);
double calculate_fitness_score(string str, string pattern);
double calculate_total_fitness_score(vector<PATTERN> best_matches);
double calculate_max_tfs(vector<PATTERN> patterns);
double calculate_min_tfs(vector<PATTERN> patterns);
int calculate_max_match(vector<PATTERN> patterns);
char select_base_pair_from_weight_matrix(WEIGHT_MATRIX weight_matrix, int index);
char select_base_pair_with_max_value(WEIGHT_MATRIX weight_matrix, int index);
char select_base_pair_without_max_value(WEIGHT_MATRIX weight_matrix, int index);
char select_base_pair_for_rearrangment(WEIGHT_MATRIX weight_matrix, int index, char present_base_pair);
pair<string, string> cross_over(string parent1, string parent2, int incision_point = 0);
pair<string, string> mutation(WEIGHT_MATRIX weight_matrix);
string rearrange_pattern(WEIGHT_MATRIX weight_matrix, string running_pattern);
string preprocess_pattern(WEIGHT_MATRIX weight_matrix);
PATTERN best_match_in_a_sequence(string sequence, string pattern);
PATTERN generate_random_pattern();
pair<vector<PATTERN>,map<string,int>> preprocess_generation1(GENERATION_INFO generation1);
WEIGHT_MATRIX generate_weight_matrix(vector<PATTERN> all_matches);
GENERATION_INFO discard_weak_patterns(GENERATION_INFO present_gener, double discarsion_factor = 0);

vector<string> sequences;
int pattern_length, number_of_pattern_in_generation1;
vector<PATTERN> best_patterns;
map<string, int> unique_best;

bool comptfs(PATTERN p1, PATTERN p2)
{
	if (p1.complete_match != p2.complete_match)
	{
		return p1.complete_match > p2.complete_match;
	}
	else
	{
		return p1.tfs > p2.tfs;
	}
	
}
bool comppat(PATTERN p1, PATTERN p2)
{
	return p1.pat > p2.pat;
}

bool check_ambiguity(string pattern)
{
	for (int i = 0; i < pattern.sz; i++)
	{
		if (pattern[i] != 'A' && pattern[i] != 'T' && pattern[i] != 'G' && pattern[i] != 'C')
		{
			return false;
		}
	}

	return true;
}

int main()
{
	ios_base::sync_with_stdio(0);
	cin.tie(NULL);

	GENERATION_INFO present_gener,previous_gener;
	char sequence[20000];
	int number_sequence, M, present_gener_no, iteration_no;

	cout << "Number of Promoter Sequence: ";
	cin >> number_sequence;

	FILE *myfile;

	myfile=fopen("test.txt","r");

	for (int i = 1; i <= number_sequence; i++)
	{
		fscanf(myfile, "%s", sequence);
		sequences.psb(sequence);
	}

	cout << "Pattern Length: ";
	cin >> pattern_length;

	cout << "Number of Pattern in Generation 1: ";
	cin >> number_of_pattern_in_generation1;

	for (int i = 1; i <= number_of_pattern_in_generation1; i++)
	{
		PATTERN pattern;

		pattern = generate_random_pattern();

		map<string, int>::iterator it;

		it = present_gener.unique_pattern.find(pattern.pat);

		if (it == present_gener.unique_pattern.end)
		{
			present_gener.unique_pattern.insert(make_pair(pattern.pat, 1));
			present_gener.patterns.psb(pattern);
		}
		else
		{
			it->second++;
		}
	}

	present_gener_no = 1;
	iteration_no = 1;

	pair<vector<PATTERN>, map<string, int>> temp;

	
	temp = preprocess_generation1(present_gener);
	present_gener.patterns = temp.fst;
	present_gener.unique_pattern = temp.snd;

	int count_same_gen = 0;
	int count_max_tfs = 0;

	while (present_gener_no <= 50)
	{
		GENERATION_INFO next_gener;

		present_gener.number_of_pattern = present_gener.patterns.sz;

		for (int i = 0; i < present_gener.number_of_pattern; i++)
		{
			PATTERN temp_pattern;
			string running_pat = present_gener.patterns[i].pat;
			present_gener.patterns[i].complete_match = 0;

			for (int j = 0; j < sequences.sz; j++)
			{
				temp_pattern = best_match_in_a_sequence(sequences[j], running_pat);
				present_gener.best_matches[running_pat].psb(temp_pattern);

				if(temp_pattern.tfs == 1.0)
				{
					present_gener.patterns[i].complete_match++;
				}
			}

			present_gener.patterns[i].tfs = calculate_total_fitness_score(present_gener.best_matches[running_pat]);
			present_gener.weight_matrix_all[running_pat] = generate_weight_matrix(present_gener.best_matches[running_pat]);
			present_gener.patterns[i].percentage=( (present_gener.patterns[i].complete_match/(sequences.sz*1.0)) * 100 );
		}

		cout << "Patterns of Generation " << present_gener_no << ": "; nl; nl;
		print_patterns(present_gener.patterns);
		nl;

		if (count_same_gen == 5 || count_max_tfs == 4)
		{
			for (int i = 0; i < present_gener.number_of_pattern; i++)
			{
				string running_pat = present_gener.patterns[i].pat;
				PATTERN temp;

				temp.pat = rearrange_pattern(present_gener.weight_matrix_all[running_pat], running_pat);
				temp.tfs = -1;

				map<string, int>::iterator it;

				it = next_gener.unique_pattern.find(temp.pat);

				if (it == next_gener.unique_pattern.end)
				{
					next_gener.unique_pattern.insert(make_pair(temp.pat, 1));
					next_gener.patterns.psb(temp);
				}
				else
				{
					it->second++;
				}
			}

			cout << "Patterns are Rearranged"; nl;
			iteration_no++;
			count_same_gen = 0;
			count_max_tfs = 0;
		}
		else
		{
			present_gener.max_tfs = calculate_max_tfs(present_gener.patterns);
			present_gener.min_tfs = calculate_min_tfs(present_gener.patterns);
			present_gener.max_match = calculate_max_match(present_gener.patterns);

			vector<PATTERN> temp = present_gener.patterns;

			sort(temp.bgn, temp.end, comptfs);

			int cnt = 0;

			for (int i = 0; i < temp.sz; i++)
			{
				if (temp[i].complete_match == 0)
				{
					break;
				}

				map<string, int>::iterator it;

				it = unique_best.find(temp[i].pat);

				if (it == unique_best.end)
				{
					unique_best.insert(make_pair(temp[i].pat, 1));
					best_patterns.psb(temp[i]);
					cnt++;
				}
				else
				{
					it->second++;
				}
			}

			present_gener = discard_weak_patterns(present_gener);

			present_gener.min_tfs = calculate_min_tfs(present_gener.patterns);
			present_gener.number_of_pattern = present_gener.patterns.sz;

			for (int i = 0; i < present_gener.number_of_pattern; i++)
			{
				PATTERN temp = present_gener.patterns[i];

				if (temp.complete_match == present_gener.max_match && temp.complete_match > 0)

				{
					map<string, int>::iterator it;

					it = next_gener.unique_pattern.find(temp.pat);

					if (it == next_gener.unique_pattern.end)
					{
						next_gener.unique_pattern.insert(make_pair(temp.pat, 1));
						next_gener.patterns.psb(temp);
					}
					else
					{
						it->second++;
					}
				}
			}
			
			for (int i = 0; i < present_gener.patterns.sz; i++)
			{
				if (present_gener.patterns[i].complete_match < present_gener.max_match || present_gener.patterns[i].complete_match == 0)
				{
					string running_pat = present_gener.patterns[i].pat;

					present_gener.parents[running_pat] = mutation(present_gener.weight_matrix_all[running_pat]);
				}
			}

			for (int i = 0; i < present_gener.patterns.sz; i++)
			{
				if (present_gener.patterns[i].complete_match < present_gener.max_match || present_gener.patterns[i].complete_match == 0)
				{
					string running_pat = present_gener.patterns[i].pat;

					present_gener.children[running_pat] = cross_over(present_gener.parents[running_pat].fst, present_gener.parents[running_pat].snd);
				}
			}

			for (int i = 0; i < present_gener.patterns.sz; i++)
			{
				if (present_gener.patterns[i].complete_match < present_gener.max_match || present_gener.patterns[i].complete_match == 0)
				{
					string running_pat = present_gener.patterns[i].pat;
					PATTERN temp;

					temp.pat = present_gener.children[running_pat].fst;
					temp.tfs = 0;

					map<string, int>::iterator it;

					it = next_gener.unique_pattern.find(temp.pat);

					if (it == next_gener.unique_pattern.end)
					{
						next_gener.unique_pattern.insert(make_pair(temp.pat, 1));
						next_gener.patterns.psb(temp);
					}
					else
					{
						it->second++;
					}

					temp.pat = present_gener.children[running_pat].snd;
					temp.tfs = 0;

					it = next_gener.unique_pattern.find(temp.pat);

					if (it == next_gener.unique_pattern.end)
					{
						next_gener.unique_pattern.insert(make_pair(temp.pat, 1));
						next_gener.patterns.psb(temp);
					}
					else
					{
						it->second++;
					}
				}
			}

			if (present_gener.unique_pattern == next_gener.unique_pattern)
			{
				count_same_gen++;
			}
			if (present_gener.max_tfs == previous_gener.max_tfs)
			{
				count_max_tfs++;
			}

			present_gener_no++;

		}

		previous_gener = present_gener;
		present_gener = next_gener;
		
	}

	print_patterns(best_patterns);

	return 0;
}

void print_patterns(vector<PATTERN> test)
{
	sort(test.bgn, test.end, comptfs);

	cout << setw(21) << "Pattern" << setw(15) << "TFS" << setw(25) << "Complete Match" << setw(20) << "Percentage"; nl; nl;

	for (int i = 0; i < test.sz; i++)
	{
		cout << setw(21) << test[i].pat << setw(16) << test[i].tfs << setw(15) << test[i].complete_match << " / " << sequences.sz << setw(23) << test[i].percentage; nl;
	}
}

double if_not_match(char base1, char base2)
{
	if ((base1 == base_pair[0] && base2 == base_pair[2]) || (base1 == base_pair[2] && base2 == base_pair[0]))
	{
		return 0.5;
	}
	else if ((base1 == base_pair[1] && base2 == base_pair[3]) || (base1 == base_pair[3] && base2 == base_pair[1]))
	{
		return 0.5;
	}
	else if ((base1 == base_pair[0] && base2 == base_pair[3]) || (base1 == base_pair[3] && base2 == base_pair[0]))
	{
		return 0.2;
	}
	else if ((base1 == base_pair[0] && base2 == base_pair[1]) || (base1 == base_pair[1] && base2 == base_pair[0]))
	{
		return 0.2;
	}
	else if ((base1 == base_pair[2] && base2 == base_pair[3]) || (base1 == base_pair[3] && base2 == base_pair[2]))
	{
		return 0.2;
	}
	else if ((base1 == base_pair[1] && base2 == base_pair[2]) || (base1 == base_pair[2] && base2 == base_pair[1]))
	{
		return 0.2;
	}
	else if (base2 == 'M' && (base1 == base_pair[0] || base1 == base_pair[3]))
	{
		return 0.2;
	}
	else if (base2 == 'R' && (base1 == base_pair[0] || base1 == base_pair[2]))
	{
		return 0.5;
	}
	else if (base2 == 'W' && (base1 == base_pair[0] || base1 == base_pair[1]))
	{
		return 0.2;
	}
	else if (base2 == 'S' && (base1 == base_pair[2] || base1 == base_pair[3]))
	{
		return 0.2;
	}
	else if (base2 == 'Y' && (base1 == base_pair[1] || base1 == base_pair[3]))
	{
		return 0.5;
	}
	else if (base2 == 'K' && (base1 == base_pair[1] || base1 == base_pair[2]))
	{
		return 0.2;
	}
	else if (base2 == 'V' && (base1 == base_pair[0] || base1 == base_pair[2] || base1 == base_pair[3]))
	{
		return 0.1;
	}
	else if (base2 == 'H' && (base1 == base_pair[0] || base1 == base_pair[1] || base1 == base_pair[3]))
	{
		return 0.1;
	}
	else if (base2 == 'D' && (base1 == base_pair[0] || base1 == base_pair[1] || base1 == base_pair[2]))
	{
		return 0.1;
	}
	else if (base2 == 'B' && (base1 == base_pair[1] || base1 == base_pair[2] || base1 == base_pair[3]))
	{
		return 0.1;
	}
	else if (base2 == 'N' && (base1 == base_pair[0] || base1 == base_pair[1] || base1 == base_pair[2] || base1 == base_pair[3]))
	{
		return 0.0;
	}
	else
	{
		return 0.0;
	}
}

double calculate_fitness_score(string str, string pattern)
{
	double result = 0;

	for (int i = 0; i < pattern.sz; i++)
	{
		if (str[i] == pattern[i])
		{
			result += 1;
		}
		else
		{
			result += if_not_match(str[i], pattern[i]);
		}
	}

	return result/pattern.sz;
}

double calculate_total_fitness_score(vector<PATTERN> best_matches)
{
	double result = 0;

	for (int i = 0; i < best_matches.sz; i++)
	{
		result += best_matches[i].tfs;
	}

	return result;
}

double calculate_max_tfs(vector<PATTERN> patterns)
{
	double temp = -1.0;

	for (int i = 0; i < patterns.sz; i++)
	{
		if (patterns[i].tfs > temp)
		{
			temp = patterns[i].tfs;
		}
	}

	return temp;
}

int calculate_max_match(vector<PATTERN> patterns)
{
	int temp = -1.0;

	for (int i = 0; i < patterns.sz; i++)
	{
		if (patterns[i].complete_match > temp)
		{
			temp = patterns[i].complete_match;
		}
	}

	return temp;
}

double calculate_min_tfs(vector<PATTERN> patterns)
{
	double temp = 1000000.0;

	for (int i = 0; i < patterns.sz; i++)
	{
		if (patterns[i].tfs < temp)
		{
			temp = patterns[i].tfs;
		}
	}

	return temp;
}

char select_base_pair_from_weight_matrix(WEIGHT_MATRIX weight_matrix, int index)
{
	if (weight_matrix.wm[0][index] > 0 && weight_matrix.wm[1][index] == 0 && weight_matrix.wm[2][index] == 0 && weight_matrix.wm[3][index] == 0)
	{
		return 'A';
	}
	else if (weight_matrix.wm[0][index] == 0 && weight_matrix.wm[1][index] > 0 && weight_matrix.wm[2][index] == 0 && weight_matrix.wm[3][index] == 0)
	{
		return 'T';
	}
	else if (weight_matrix.wm[0][index] == 0 && weight_matrix.wm[1][index] == 0 && weight_matrix.wm[2][index] > 0 && weight_matrix.wm[3][index] == 0)
	{
		return 'G';
	}
	else if (weight_matrix.wm[0][index] == 0 && weight_matrix.wm[1][index] == 0 && weight_matrix.wm[2][index] == 0 && weight_matrix.wm[3][index] > 0)
	{
		return 'C';
	}
	else if (weight_matrix.wm[0][index] > 0 && weight_matrix.wm[1][index] > 0 && weight_matrix.wm[2][index] > 0 && weight_matrix.wm[3][index] > 0)
	{
		return 'N';
	}
	else if (weight_matrix.wm[0][index] > 0 && weight_matrix.wm[1][index] > 0 && weight_matrix.wm[2][index] > 0 && weight_matrix.wm[3][index] == 0)
	{
		return 'D';
	}
	else if (weight_matrix.wm[0][index] > 0 && weight_matrix.wm[1][index] > 0 && weight_matrix.wm[2][index] == 0 && weight_matrix.wm[3][index] > 0)
	{
		return 'H';
	}
	else if (weight_matrix.wm[0][index] > 0 && weight_matrix.wm[1][index] == 0 && weight_matrix.wm[2][index] > 0 && weight_matrix.wm[3][index] > 0)
	{
		return 'V';
	}
	else if (weight_matrix.wm[0][index] == 0 && weight_matrix.wm[1][index] > 0 && weight_matrix.wm[2][index] > 0 && weight_matrix.wm[3][index] > 0)
	{
		return 'B';
	}
	else if (weight_matrix.wm[0][index] > 0 && weight_matrix.wm[1][index] == 0 && weight_matrix.wm[2][index] == 0 && weight_matrix.wm[3][index] > 0)
	{
		return 'M';
	}
	else if (weight_matrix.wm[0][index] > 0 && weight_matrix.wm[1][index] == 0 && weight_matrix.wm[2][index] > 0 && weight_matrix.wm[3][index] == 0)
	{
		return 'R';
	}
	else if (weight_matrix.wm[0][index] > 0 && weight_matrix.wm[1][index] > 0 && weight_matrix.wm[2][index] == 0 && weight_matrix.wm[3][index] == 0)
	{
		return 'W';
	}
	else if (weight_matrix.wm[0][index] == 0 && weight_matrix.wm[1][index] == 0 && weight_matrix.wm[2][index] > 0 && weight_matrix.wm[3][index] > 0)
	{
		return 'S';
	}
	else if (weight_matrix.wm[0][index] == 0 && weight_matrix.wm[1][index] > 0 && weight_matrix.wm[2][index] == 0 && weight_matrix.wm[3][index] > 0)
	{
		return 'Y';
	}
	else if (weight_matrix.wm[0][index] == 0 && weight_matrix.wm[1][index] > 0 && weight_matrix.wm[2][index] > 0 && weight_matrix.wm[3][index] == 0)
	{
		return 'K';
	}
}

char select_base_pair_with_max_value(WEIGHT_MATRIX weight_matrix, int index)
{
	double mx = max4(weight_matrix.wm[0][index], weight_matrix.wm[1][index], weight_matrix.wm[2][index], weight_matrix.wm[3][index]);

	for (int j = 0; j < 4; j++)
	{
		if (weight_matrix.wm[j][index] != mx)
		{
			weight_matrix.wm[j][index] = 0;
		}
	}

	return select_base_pair_from_weight_matrix(weight_matrix, index);
}

char select_base_pair_without_max_value(WEIGHT_MATRIX weight_matrix, int index)
{
	double mx = max4(weight_matrix.wm[0][index], weight_matrix.wm[1][index], weight_matrix.wm[2][index], weight_matrix.wm[3][index]);

	for (int j = 0; j < 4; j++)
	{
		if (weight_matrix.wm[j][index] == mx && mx != 1)
		{
			weight_matrix.wm[j][index] = 0;
			break;
		}
	}

	return select_base_pair_with_max_value(weight_matrix, index);
}

char select_base_pair_for_rearrangment(WEIGHT_MATRIX weight_matrix, int index, char present_base_pair)
{
	double mx = max4(weight_matrix.wm[0][index], weight_matrix.wm[1][index], weight_matrix.wm[2][index], weight_matrix.wm[3][index]);

	for (int j = 0; j < 4; j++)
	{
		if (mx != 1 && base_pair[j] == present_base_pair)
		{
			weight_matrix.wm[j][index] = 0;
			break;
		}
	}

	return select_base_pair_with_max_value(weight_matrix, index);
}

pair<string, string> cross_over(string parent1, string parent2, int incision_point)
{
	if (incision_point == 0)
	{
		incision_point = parent1.sz / 2;
	}

	string p1l, p1r, p2l, p2r, child1, child2;

	p1l.assign(parent1.bgn, parent1.bgn + incision_point);
	p1r.assign(parent1.bgn + incision_point, parent1.end);
	p2l.assign(parent2.bgn, parent2.bgn + incision_point);
	p2r.assign(parent2.bgn + incision_point, parent2.end);

	child1 = (p1l + p2r);
	child2 = (p2l + p1r);

	return make_pair(child1, child2);
}

pair<string, string> mutation(WEIGHT_MATRIX weight_matrix)
{
	string parent1, parent2;

	for (int i = 0; i < weight_matrix.wm[0].sz; i++)
	{
		parent1.psb(select_base_pair_with_max_value(weight_matrix, i));
		parent2.psb(select_base_pair_without_max_value(weight_matrix, i));
	}

	return make_pair(parent1, parent2);
}

string rearrange_pattern(WEIGHT_MATRIX weight_matrix,string running_pattern)
{
	string rearranged_pattern;

	for (int i = 0; i < pattern_length; i++)
	{
		rearranged_pattern.psb(select_base_pair_for_rearrangment(weight_matrix, i,running_pattern[i]));
	}

	return rearranged_pattern;
}

string preprocess_pattern(WEIGHT_MATRIX weight_matrix)
{
	string preprocessed_pattern;

	for (int i = 0; i < weight_matrix.wm[0].sz; i++)
	{
		preprocessed_pattern.psb(select_base_pair_from_weight_matrix(weight_matrix, i));
	}

	return preprocessed_pattern;
}

PATTERN best_match_in_a_sequence(string sequence, string pattern)
{
	double mx;
	PATTERN best_pattern;
	best_pattern.tfs = -1.0;

	for (int i = 0; i < sequence.sz - pattern.sz; i++)
	{
		string temp(sequence.bgn + i, sequence.bgn + i + pattern.sz);

		mx = calculate_fitness_score(temp, pattern);

		if (mx > best_pattern.tfs)
		{
			best_pattern.tfs = mx;
			best_pattern.pat = temp;
		}
	}

	return best_pattern;
}

PATTERN generate_random_pattern()
{
	PATTERN generated_pattern;

	for (int i = 0; i < pattern_length; i++)
	{
		int x = rand() % sequences.sz;
		int y = rand() % sequences[x].sz;

		generated_pattern.pat.psb(sequences[x][y]);
	}

	generated_pattern.tfs = -1;

	return generated_pattern;
}

pair<vector<PATTERN>,map<string,int>> preprocess_generation1(GENERATION_INFO generation1)
{
	map<string, int> temp;
	PATTERN temppat;
	vector<PATTERN> preprocessed;

	for (int i = 0; i < generation1.patterns.sz; i++)
	{
		PATTERN temp_pattern;

		for (int j = 0; j < sequences.sz; j++)
		{
			temp_pattern = best_match_in_a_sequence(sequences[j], generation1.patterns[i].pat);
			generation1.best_matches[generation1.patterns[i].pat].psb(temp_pattern);
		}

		generation1.weight_matrix_all[generation1.patterns[i].pat] = generate_weight_matrix(generation1.best_matches[generation1.patterns[i].pat]);
		
		temppat.pat = preprocess_pattern(generation1.weight_matrix_all[generation1.patterns[i].pat]);

		map<string, int>::iterator it;

		it = temp.find(temppat.pat);

		if (it == temp.end)
		{
			temppat.tfs = -1;
			preprocessed.psb(temppat);
			temp.insert(make_pair(temppat.pat, 1));
		}
		else
		{
			it->snd++;
		}
	}

	return make_pair(preprocessed,temp);
}

WEIGHT_MATRIX generate_weight_matrix(vector<PATTERN> all_matches)
{
	WEIGHT_MATRIX weight_matrix;
	int cnt_t, cnt_g, cnt_c, cnt_a;

	for (int i = 0; i < pattern_length; i++)
	{
		cnt_a = cnt_c = cnt_g = cnt_t = 0;

		for (int j = 0; j < all_matches.sz; j++)
		{
			if (all_matches[j].pat[i] == base_pair[0])
			{
				cnt_a++;
			}
			else if (all_matches[j].pat[i] == base_pair[1])
			{
				cnt_t++;
			}
			else if (all_matches[j].pat[i] == base_pair[2])
			{
				cnt_g++;
			}
			else if (all_matches[j].pat[i] == base_pair[3])
			{
				cnt_c++;
			}
		}

		weight_matrix.wm[0].psb((cnt_a*1.0) / (all_matches.sz*1.0));
		weight_matrix.wm[1].psb((cnt_t*1.0) / (all_matches.sz*1.0));
		weight_matrix.wm[2].psb((cnt_g*1.0) / (all_matches.sz*1.0));
		weight_matrix.wm[3].psb((cnt_c*1.0) / (all_matches.sz*1.0));

	}

	return weight_matrix;
}

GENERATION_INFO discard_weak_patterns(GENERATION_INFO present_gener, double discarsion_factor)
{
	if (present_gener.number_of_pattern <=  (2*number_of_pattern_in_generation1))
	{
		discarsion_factor = (present_gener.max_tfs - present_gener.min_tfs) / 2.5;
	}
	else
	{
		discarsion_factor = (present_gener.max_tfs - present_gener.min_tfs) / 7;
	}

	int cnt = 0;

	for (int i = 0; i < present_gener.patterns.sz;)
	{
		PATTERN running_pat = present_gener.patterns[i];

		if ((running_pat.tfs < present_gener.max_tfs - discarsion_factor && check_ambiguity(running_pat.pat)) || cnt >= 2*number_of_pattern_in_generation1)
		{
			map<string, vector<PATTERN> >::iterator it1;
			map<string, int>::iterator it2;
			map<string, WEIGHT_MATRIX>::iterator it3;

			it1 = present_gener.best_matches.find(running_pat.pat);

			if (it1 != present_gener.best_matches.end)
			{
				present_gener.best_matches.erase(it1);

			}
			it2 = present_gener.unique_pattern.find(running_pat.pat);

			if (it2 != present_gener.unique_pattern.end)
			{
				present_gener.unique_pattern.erase(it2);

			}

			it3 = present_gener.weight_matrix_all.find(running_pat.pat);

			if (it3 != present_gener.weight_matrix_all.end)
			{
				present_gener.weight_matrix_all.erase(it3);
			}

			present_gener.patterns.erase(present_gener.patterns.bgn + i);
		}
		else
		{
			i++;
			cnt++;
		}
	}

	return present_gener;
}
