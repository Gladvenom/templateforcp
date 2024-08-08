//if (strcmp(mpp[pattern[i]].c_str(), curw.c_str()) != 0) for comparing string
string matching algo (rabin karp)
//  vector<int> rabin_karp(string const& s, string const& t) {
//         const int p = 31; 
//         const int m = 1e9 + 21;
//         int S = s.size(), T = t.size();

//         vector<long long> p_pow(max(S, T)); 
//         p_pow[0] = 1; 
//         for (int i = 1; i < (int)p_pow.size(); i++) 
//             p_pow[i] = (p_pow[i-1] * p) % m;

//         vector<long long> h(T + 1, 0); 
//         for (int i = 0; i < T; i++)
//             h[i+1] = (h[i] + (t[i] - 'a' + 1) * p_pow[i]) % m; 
//         long long h_s = 0; 
//         for (int i = 0; i < S; i++) 
//             h_s = (h_s + (s[i] - 'a' + 1) * p_pow[i]) % m; 

//         vector<int> occurrences;
//         for (int i = 0; i + S - 1 < T; i++) {
//             long long cur_h = (h[i+S] + m - h[i]) % m;
//             if (cur_h == h_s * p_pow[i] % m)
//                 occurrences.push_back(i);
//         }
//         return occurrences;
//     }
 string matching algo (KMP)
// std::vector<int> computeLPS(const std::string& pattern) {
//     int m = pattern.length();
//     std::vector<int> lps(m, 0);
//     int len = 0;  // Length of the previous longest prefix suffix

//     for (int i = 1; i < m; ) {
//         if (pattern[i] == pattern[len]) {
//             len++;
//             lps[i] = len;
//             i++;
//         } else {
//             if (len != 0) {
//                 len = lps[len - 1];
//             } else {
//                 lps[i] = 0;
//                 i++;
//             }
//         }
//     }
//     return lps;
// }

// // Function to find occurrences of pattern in text using KMP algorithm
// std::vector<int> kmpSearch(const std::string& text, const std::string& pattern) {
//     int n = text.length();
//     int m = pattern.length();
//     std::vector<int> lps = computeLPS(pattern);
//     std::vector<int> positions;

//     int i = 0;  // Index for text[]
//     int j = 0;  // Index for pattern[]

//     while (i < n) {
//         if (pattern[j] == text[i]) {
//             j++;
//             i++;
//         }
//         if (j == m) {
//             positions.push_back(i - j);
//             j = lps[j - 1];
//         } else if (i < n && pattern[j] != text[i]) {
//             if (j != 0)
//                 j = lps[j - 1];
//             else
//                 i++;
//         }
//     }
//     return positions;
// }
String hash(for matching string using map)
 #define ll long long
    //static vector<ll>kk,primes;
 
// change '0' to '#' if needed or remove '0';
// struct string_hash {
//   int len;
//   ll mod, poly, inv;
//   vector<ll> prefix;
//   vector<ll> invs;

//   void init(int n, string s, ll k = 89, ll m = 1000000007) {
//     mod = m;
//     poly = k;
//     prefix = vector<ll>(n);
//     invs = vector<ll>(n);

//     invs[0] = 1;
//     inv = minv(k);
//     for (int i = 1; i < n; i++) {
//       invs[i] = (invs[i - 1] * inv) % mod;
//     }

//     ll x = 1;
//     prefix[0] = (s[0] - '0' + 1);
//     for (int i = 1; i < n; i++) {
//       x = (x * k) % mod;
//       prefix[i] = (prefix[i - 1] + x * (s[i] - '0' + 1)) % mod;
//     }

//     len = n;
//   }

//   void extend(string next) {
//     int x = next.length();
//     for (int i = 0; i < x; i++) {
//       invs.push_back((invs[i - 1] * inv) % mod);
//     }

//     ll p = mpow(poly, len - 1);
//     for (int i = 0; i < x; i++) {
//       p = (p * poly) % mod;
//       prefix.push_back((prefix[i - 1] + p * (next[i - len] - '0' + 1)) % mod);
//     }

//     len += x;
//   }

//   ll get_hash(int left, int right) {
//     if (left == 0) return prefix[right];
//     return ((prefix[right] - prefix[left - 1] + mod) * invs[left]) % mod;
//   }

//   ll mpow(ll base, ll exp) {
//     ll res = 1;
//     while (exp) {
//       if (exp % 2 == 1) {
//         res = (res * base) % mod;
//       }
//       exp >>= 1;
//       base = (base * base) % mod;
//     }
//     return res;
//   }
//   ll minv(ll base) {
//     return mpow(base, mod - 2);
//   }
// };

// template<int K>
// struct multihash {
//   string_hash sh[K];
//     vector<ll>kk {
//   89,
//   101,
//   189,
//   94,
//   621,
//   (ll)rand() % 1000 + 101,
//   (ll)rand() % 2000 + 121,
//   (ll)rand() % 4000 + 141,
//   (ll)rand() % 8000 + 161,
//   (ll)rand() % 16000 + 183
// };

//  vector<ll>primes {
//   533000401,
//   735632791,
//   776531419,
//   797003413,
//   817504243,
//   920419813,
//   961748941,
//   982451653,
//   1000000007,
//   1000000009
// };

//   void init(int n, string s) {
       
//     for (int i = 0; i < K; i++) {
//       sh[i].init(n, s, kk[i], primes[9 - i]);
//     }
//   }

//   vector<ll> get_hash(int l, int r) {
//     vector<ll> ret(K);
//     for (int i = 0; i < K; i++) {
//       ret[i] = sh[i].get_hash(l, r);
//     }
//     return ret;
//   }

//   bool compare(vector<ll> a, vector<ll> b) {
//     for (int i = 0; i < K; i++) {
//       if (a[i] != b[i]) return 0;
//     }
//     return 1;
//   }
// };
string matching(Z function)
//     vector<int> z_function(string s) {
//     int n = (int) s.length();
//     vector<int> z(n);
//     for (int i = 1, l = 0, r = 0; i < n; ++i) {
//         if (i <= r)
//             z[i] = min (r - i + 1, z[i - l]);
//         while (i + z[i] < n && s[z[i]] == s[i + z[i]])
//             ++z[i];
//         if (i + z[i] - 1 > r)
//             l = i, r = i + z[i] - 1;
//     }
//     z[0]=s.size();
//     return z;
// }
segment tree template striver
// class SGTree {
// 	vector<int> seg;
// public:
// 	SGTree(int n) {
// 		seg.resize(4 * n + 1);
// 	}

// 	void build(int ind, int low, int high, int arr[]) {
// 		if (low == high) {
// 			seg[ind] = arr[low];
// 			return;
// 		}

// 		int mid = (low + high) / 2;
// 		build(2 * ind + 1, low, mid, arr);
// 		build(2 * ind + 2, mid + 1, high, arr);
// 		seg[ind] = min(seg[2 * ind + 1], seg[2 * ind + 2]);
// 	}

// 	int query(int ind, int low, int high, int l, int r) {
// 		// no overlap
// 		// l r low high or low high l r
// 		if (r < low || high < l) return INT_MAX;

// 		// complete overlap
// 		// [l low high r]
// 		if (low >= l && high <= r) return seg[ind];

// 		int mid = (low + high) >> 1;
// 		int left = query(2 * ind + 1, low, mid, l, r);
// 		int right = query(2 * ind + 2, mid + 1, high, l, r);
// 		return min(left, right);
// 	}
// 	void update(int ind, int low, int high, int i, int val) {
// 		if (low == high) {
// 			seg[ind] = val;
// 			return;
// 		}

// 		int mid = (low + high) >> 1;
// 		if (i <= mid) update(2 * ind + 1, low, mid, i, val);
// 		else update(2 * ind + 2, mid + 1, high, i, val);
// 		seg[ind] = min(seg[2 * ind + 1], seg[2 * ind + 2]);
// 	}
// };
N C r in O(1)
// #include <bits/stdc++.h> 
// #define int  long long 

// const int N = 1000001; 

// using namespace std; 
  
// // array to store inverse of 1 to N 
// int  factorialNumInverse[N + 1]; 
  
// // array to precompute inverse of 1! to N! 
// int  naturalNumInverse[N + 1]; 
  
// // array to store factorial of first N numbers 
// int  fact[N + 1]; 
  
// // Function to precompute inverse of numbers 
// void InverseofNumber(int  p) 
// { 
//     naturalNumInverse[0] = naturalNumInverse[1] = 1; 
//     for (int i = 2; i <= N; i++) 
//         naturalNumInverse[i] = naturalNumInverse[p % i] * (p - p / i) % p; 
// } 
// // Function to precompute inverse of factorials 
// void InverseofFactorial(int  p) 
// { 
//     factorialNumInverse[0] = factorialNumInverse[1] = 1; 
  
//     // precompute inverse of natural numbers 
//     for (int i = 2; i <= N; i++) 
//         factorialNumInverse[i] = (naturalNumInverse[i] * factorialNumInverse[i - 1]) % p; 
// } 
  
// // Function to calculate factorial of 1 to N 
// void factorial(int  p) 
// { 
//     fact[0] = 1; 
  
//     // precompute factorials 
//     for (int i = 1; i <= N; i++) { 
//         fact[i] = (fact[i - 1] * i) % p; 
//     } 
// } 
  
// // Function to return nCr % p in O(1) time 
// int  Binomial(int  N, int  R, int  p) 
// { 
//     // n C r = n!*inverse(r!)*inverse((n-r)!) 
//     int  ans = ((fact[N] * factorialNumInverse[R]) 
//               % p * factorialNumInverse[N - R]) 
//              % p; 
//     return ans; 
// } 
  

// int main() 
// { 
   
//     int  p = 1000000007; 
//     InverseofNumber(p); 
//     InverseofFactorial(p); 
//     factorial(p); 
  
//     // 1st query 
//     int  N = 15; 
//     int  R = 4; 
//     cout << Binomial(N, R, p) << endl; 
  
//     // 2nd query 
//     N = 20; 
//     R = 3; 
//     cout << Binomial(N, R, p) << endl; 
  
//     return 0; 
// } 
