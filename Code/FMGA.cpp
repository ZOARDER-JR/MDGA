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
char base_pair_amb[] = { 'M', 'R', 'W', 'S', 'Y', 'K' };

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
	int number_of_pattern;
	double max_tfs, min_tfs;
	vector<PATTERN> patterns;
	vector<vector<PATTERN>> best_matches;
	vector<WEIGHT_MATRIX> weight_matrix_all;
	vector<pair<string, string>> parents;
	vector<pair<string, string>> children;
};

void print_patterns(vector<PATTERN> test);
double if_not_match(char base1, char base2);
double calculate_fitness_score(string str, string pattern);
double calculate_total_fitness_score(vector<PATTERN> best_matches);
double calculate_max_tfs(vector<PATTERN> patterns);
double calculate_min_tfs(vector<PATTERN> patterns);
double calcutele_ambiguity_code_penalty(string pattern);
char select_base_pair_from_weight_matrix(WEIGHT_MATRIX weight_matrix, int index);
char select_base_pair_with_max_value(WEIGHT_MATRIX weight_matrix, int index);
char select_base_pair_without_max_value(WEIGHT_MATRIX weight_matrix, int index);
char select_base_pair_for_rearrangment(WEIGHT_MATRIX weight_matrix, int index, char present_base_pair);
pair<string, string> mutation(WEIGHT_MATRIX weight_matrix);
pair<string, string> cross_over(string parent1, string parent2, int incision_point=0);
string rearrange_pattern(WEIGHT_MATRIX weight_matrix);
PATTERN best_match_in_a_sequence(string sequence, string pattern);
PATTERN generate_random_pattern();
WEIGHT_MATRIX generate_weight_matrix(vector<PATTERN> all_matches);

vector<string> sequences;
int pattern_length;
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
bool isSame(vector<PATTERN> v1, vector<PATTERN> v2)
{
	sort(v1.bgn, v1.end, comppat);
	sort(v2.bgn, v2.end, comppat);

	for (int i = 0; i < v1.sz; i++)
	{
		if (v1[i].pat != v2[i].pat)
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

	GENERATION_INFO present_gener, previous_gener;
	char sequence[20000];
	int number_sequence, number_of_pattern_in_generation1, M, present_gener_no, iteration_no;

	cout << "Number of Promoter Sequence: ";
	cin >> number_sequence;

	FILE *myfile;

	myfile = fopen("test.txt", "r");

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

		present_gener.patterns.psb(pattern);
	}

	present_gener_no = 0;

	cout << "Number of Highest Iteration: ";
	cin >> M;
	iteration_no = 1;

	int count_same_gen = 0;

	while (iteration_no <= M)
	{
		present_gener_no++;

		if (present_gener_no > 50)
		{
			break;
		}

		present_gener.number_of_pattern = present_gener.patterns.sz;

		vector<PATTERN> T;
		WEIGHT_MATRIX W;
		pair<string, string> P;

		present_gener.best_matches.assign(present_gener.number_of_pattern,T);
		present_gener.weight_matrix_all.assign(present_gener.number_of_pattern, W);
		present_gener.parents.assign(present_gener.number_of_pattern, P);
		present_gener.children.assign(present_gener.number_of_pattern, P);

		GENERATION_INFO next_gener;

		for (int i = 0; i < present_gener.number_of_pattern; i++)
		{
			PATTERN temp_pattern;
			string running_pat = present_gener.patterns[i].pat;
			present_gener.patterns[i].complete_match = 0;

			for (int j = 0; j < sequences.sz; j++)
			{
				temp_pattern = best_match_in_a_sequence(sequences[j], running_pat);
				present_gener.best_matches[i].psb(temp_pattern);
				if( temp_pattern.tfs == 1)
				{
					present_gener.patterns[i].complete_match++;
				}
			}

			present_gener.weight_matrix_all[i] = generate_weight_matrix(present_gener.best_matches[i]);
			present_gener.patterns[i].percentage=( (present_gener.patterns[i].complete_match/(sequences.sz*1.0)) * 100 );
		}

		for (int i = 0; i < present_gener.number_of_pattern; i++)
		{
			string running_pat = present_gener.patterns[i].pat;

			present_gener.patterns[i].tfs = calculate_total_fitness_score(present_gener.best_matches[i]);
		}

		present_gener.max_tfs = calculate_max_tfs(present_gener.patterns);
		present_gener.min_tfs = calculate_min_tfs(present_gener.patterns);

		if (count_same_gen == 5)
		{
			for (int i = 0; i < present_gener.number_of_pattern; i++)
			{
				string running_pat = present_gener.patterns[i].pat;
				PATTERN temp;

				temp.pat = rearrange_pattern(present_gener.weight_matrix_all[i]);

				temp.tfs = -1;

				next_gener.patterns.psb(temp);
			}

			cout << "Patterns are Rearranged"; nl;
			iteration_no++;
			count_same_gen = 0;
		}
		else
		{
			vector<PATTERN> temp = present_gener.patterns;

			sort(temp.bgn, temp.end, comptfs);

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
				}
				else
				{
					it->second++;
				}
			}

			for (int i = 0; i < present_gener.patterns.sz; i++)
			{
				if (present_gener.patterns[i].tfs == present_gener.max_tfs)
				{
					next_gener.patterns.psb(present_gener.patterns[i]);
				}
			}

			for (int i = 0; i < present_gener.number_of_pattern; i++)
			{
				if (present_gener.patterns[i].tfs < present_gener.max_tfs)
				{
					string running_pat = present_gener.patterns[i].pat;

					present_gener.parents[i] = mutation(present_gener.weight_matrix_all[i]);
				}
			}

			for (int i = 0; i < present_gener.patterns.sz; i++)
			{
				if (present_gener.patterns[i].tfs < present_gener.max_tfs)
				{
					string running_pat = present_gener.patterns[i].pat;

					present_gener.children[i] = cross_over(present_gener.parents[i].fst, present_gener.parents[i].snd);
				}
			}

			for (int i = 0; i < present_gener.patterns.sz; i++)
			{
				if (present_gener.patterns[i].tfs < present_gener.max_tfs)
				{
					string running_pat = present_gener.patterns[i].pat;
					PATTERN temp, child1, child2;
					vector<PATTERN> c1_best_matches, c2_best_matches;

					child1.pat = present_gener.children[i].fst;
					child2.pat = present_gener.children[i].snd;

					double totalc1;
					double totalc2;

					totalc1 = 0.0;
					totalc2 = 0.0;

					for (int j = 0; j < pattern_length; j++)
					{
						for (int k = 0; k < 4; k++)
						{
							if (child1.pat[j] == base_pair[k])
							{
								totalc1 += 1;
								break;
							}
							if (child2.pat[j] == base_pair[k])
							{
								totalc2 += 1;
								break;
							}
						}
						for (int k = 0; k < 6; k++)
						{
							if (child1.pat[j] == base_pair_amb[k])
							{
								totalc1 += .5;
							}
							if (child2.pat[j] == base_pair_amb[k])
							{
								totalc2 += .5;
							}
						}
					}

					child1.tfs = totalc1 - calcutele_ambiguity_code_penalty(child1.pat);
					child2.tfs = totalc2 - calcutele_ambiguity_code_penalty(child2.pat);

					temp = (child1.tfs > child2.tfs ? child1 : child2);
					temp.tfs = -1;

					next_gener.patterns.psb(temp);
				}

			}

			if (isSame(present_gener.patterns,next_gener.patterns))
			{
				count_same_gen++;
			}
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

double calcutele_ambiguity_code_penalty(string pattern)
{
	double aco, ep1, ep2, result;

	aco = ep1 = ep2 = result = 0.0;

	for (int j = 0; j < pattern_length; j++)
	{
		for (int k = 0; k < 6; k++)
		{
			if (pattern[j] == base_pair_amb[k])
			{
				ep1 += 1;

				int x = j + 1;
				while (x < pattern_length)
				{
					if (pattern[x] == 'A' || pattern[x] == 'T' || pattern[x] == 'G' || pattern[x] == 'C')
					{
						aco++;
						break;
					}

					x++;
				}
			}
		}
		if (pattern[j] == 'N')
		{
			ep2 += 1;

			int x = j + 1;
			while (x < pattern_length)
			{
				if (pattern[x] == 'A' || pattern[x] == 'T' || pattern[x] == 'G' || pattern[x] == 'C')
				{
					aco++;
					break;
				}

				x++;
			}
		}
	}

	result = ((aco * .5) + (ep1*.3) + (ep2*.5));
	return result;
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
		return 0.5;
	}
	else if ((base1 == base_pair[0] && base2 == base_pair[1]) || (base1 == base_pair[1] && base2 == base_pair[0]))
	{
		return 0.5;
	}
	else if ((base1 == base_pair[2] && base2 == base_pair[3]) || (base1 == base_pair[3] && base2 == base_pair[2]))
	{
		return 0.5;
	}
	else if ((base1 == base_pair[1] && base2 == base_pair[2]) || (base1 == base_pair[2] && base2 == base_pair[1]))
	{
		return 0.5;
	}
	else if (base2 == 'M' && (base1 == base_pair[0] || base1 == base_pair[3]))
	{
		return 0.5;
	}
	else if (base2 == 'R' && (base1 == base_pair[0] || base1 == base_pair[2]))
	{
		return 0.5;
	}
	else if (base2 == 'W' && (base1 == base_pair[0] || base1 == base_pair[1]))
	{
		return 0.5;
	}
	else if (base2 == 'S' && (base1 == base_pair[2] || base1 == base_pair[3]))
	{
		return 0.5;
	}
	else if (base2 == 'Y' && (base1 == base_pair[1] || base1 == base_pair[3]))
	{
		return 0.5;
	}
	else if (base2 == 'K' && (base1 == base_pair[1] || base1 == base_pair[2]))
	{
		return 0.5;
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

	result -= calcutele_ambiguity_code_penalty(pattern);
	return (result/pattern.sz);
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

double calculate_min_tfs(vector<PATTERN> patterns)
{
	double temp = 2.0;

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

string rearrange_pattern(WEIGHT_MATRIX weight_matrix)
{
	string rearranged_pattern;

	for (int i = 0; i < pattern_length; i++)
	{
		rearranged_pattern.psb(select_base_pair_with_max_value(weight_matrix, i));
	}

	return rearranged_pattern;
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

	generated_pattern.tfs = 0;

	return generated_pattern;
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