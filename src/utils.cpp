#include "utils.h"

/**
 * Splits a line into a vector - PERL style
 */
void split(std::vector<std::string>& vec, const char *delims, std::string& str, uint32_t limit, bool clear, bool collapse)
{
    std::map<char, int32_t> delim_set;

    for (uint32_t i=0; i<strlen(delims); ++i)
    {
        delim_set[delims[i]] = 1;
    }

    if (clear)
    {
        vec.clear();
    }
    const char* tempStr = str.c_str();
    int32_t i=0, lastIndex = str.size()-1;
    std::stringstream token;

    if (lastIndex<0) return;

    uint32_t noTokens = 0;
    bool isDelim = false;
    while (i<=lastIndex)
    {
        isDelim = (delim_set.find(tempStr[i])!=delim_set.end());

        if (!isDelim || noTokens>=limit-1)
        {
            token << tempStr[i];
        }

        if ((isDelim && noTokens<limit-1) || i==lastIndex)
        {
            if (collapse && token.str()!="")
            {
                vec.push_back(token.str());
                ++noTokens;
                token.str("");
            }
        }

        ++i;
    }
};

/**
 * Splits a line into a vector - PERL style
 */
void split(std::vector<std::string>& vec, const char *delims, const char* str, uint32_t limit, bool clear, bool collapse)
{
    std::map<char, int32_t> delim_set;

    for (uint32_t i=0; i<strlen(delims); ++i)
    {
        delim_set[delims[i]] = 1;
    }

    if (clear)
    {
        vec.clear();
    }
    const char* tempStr = str;
    int32_t i=0, lastIndex = strlen(str)-1;
    std::stringstream token;

    if (lastIndex<0) return;

    uint32_t noTokens = 0;
    bool isDelim = false;
    while (i<=lastIndex)
    {
        isDelim = (delim_set.find(tempStr[i])!=delim_set.end());

        if (!isDelim || noTokens>=limit-1)
        {
            token << tempStr[i];
        }

        if ((isDelim && noTokens<limit-1) || i==lastIndex)
        {
            if (collapse && token.str()!="")
            {
                vec.push_back(token.str());
                ++noTokens;
                token.str("");
            }
        }

        ++i;
    }
};

/**
 * Casts a string into int32.  Returns true if successful.
 */
bool str2int32(std::string& s, int32_t& i)
{
    const char* start = s.c_str();
    char *end = 0;
    i = std::strtol(s.c_str(), &end, 10);
    return (end!=start);
};

/**
 * Casts a string into uint32.  Returns true if successful.
 */
bool str2uint32(std::string& s, uint32_t& i)
{
    const char* start = s.c_str();
    char *end = 0;
    i = std::strtoul(s.c_str(), &end, 10);
    return (end!=start);
};

/**
 * Appends cuurent working directoy to a path.
 * Returns true if successful.
 */
bool append_cwd(std::string& path)
{
    if (path.size()>0 && path.c_str()[0]!='/')
    {
        char cwd[1024];
        if (getcwd(cwd, sizeof(cwd))!=NULL)
        {
            std::string cwd_path(cwd);
            path = cwd_path + "/" + path;

            return true;
        }
    }

    return false;
};


unsigned int str_hash(const char* s, unsigned int seed)
{
  unsigned int hash = seed;
  while (*s) {
    hash = hash * 101  +  *s++;
  }
  return hash;
};


/**
 * Get the last line of a file .
 */
std::string GetLastLine(const std::string& f) {
    std::ifstream file(f);
    file.seekg(-1, std::ios_base::end);
    char c;
    file.get(c);
    while (c == '\n') {
      file.seekg(-2,std::ios_base::cur);
      file.get(c);
    }
    while (c != '\n') {
      file.seekg(-2,std::ios_base::cur);
      file.get(c);
    }
    std::string line;
    std::getline(file, line);
    return line;
};


/**
 * Generate all possible k-mer from an alphabet.
 */
void EnumerateMotif(std::vector<char>& alphabet, int32_t k, std::vector<std::string >& res) {
  uint32_t n = alphabet.size();
  if (k == 1) {
    for (uint32_t i = 0; i < n; ++i) {
      std::string tmp(1,alphabet[i]);
      res.push_back(tmp);
    }
    return;
  }
  std::vector<std::string > sres;
  EnumerateMotif(alphabet, k-1, sres);
  for (uint32_t i = 0; i < n; ++i) {
    std::vector<std::string > nres = sres;
    for (uint32_t j = 0; j < nres.size(); j++) {
      nres[j] += alphabet[i];
    }
    res.insert(res.end(), nres.begin(), nres.end());
  }
  return;
};

int32_t AllConfig(std::vector<char>& alphabet, int32_t k, std::vector<std::string >& res) {
  EnumerateMotif(alphabet,k,res);
  return res.size();
};

/**
 * Sample without replacement.
 */
void NchooseK(int32_t N, int32_t k, std::set<int32_t> & chosen, std::mt19937& rng, int32_t base) {
  if (base == 0) {
    if (N == k) {
      for (int32_t r = 0; r < N; ++r) {chosen.insert(r);}
    } else {
      for (int32_t r = N - k - 1; r < N-1; ++r) {
        int32_t v = std::uniform_int_distribution<>(0, r)(rng);
        if (!chosen.insert(v).second) {
          chosen.insert(r+1);
        }
      }
    }
  } else {
    if (N == k) {
      for (int32_t r = 1; r <= N; ++r) {chosen.insert(r);}
    } else {
      for (int32_t r = N - k; r < N; ++r) {
        int32_t v = std::uniform_int_distribution<>(1, r)(rng);
        if (!chosen.insert(v).second) {
          chosen.insert(r+1);
        }
      }
    }
  }
};

void NchooseK(int32_t N, int32_t k, std::set<int32_t> & chosen, std::vector<int32_t> & avoid, std::mt19937& rng, int32_t base) {
  if (base == 0) {
    for (int32_t r = N - k - 1; r < N-((int32_t) avoid.size())-1; ++r) {
      int32_t v = std::uniform_int_distribution<>(0, r)(rng);
      if (!chosen.insert(v).second) {
        chosen.insert(r+1);
      }
    }
  } else {
    for (int32_t r = N - k; r < N; ++r) {
      int32_t v = std::uniform_int_distribution<>(1, r)(rng);
      if (!chosen.insert(v).second) {
        chosen.insert(r+1);
      }
    }
  }
  int32_t ct = 1 - base;
  for (auto & v : avoid) {
    auto it = chosen.find(v);
    if (it != chosen.end()) {
      chosen.erase(it);
      chosen.insert(N-ct);
    }
    ct++;
  }
};










