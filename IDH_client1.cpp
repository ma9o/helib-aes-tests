/** Test_AES.cpp - test program for homomorphic AES using HElib
 */
#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

namespace std {} using namespace std;
namespace NTL {} using namespace NTL;
namespace helib{} using namespace helib;
#include <cstring>
#include "homAES.h"
#include <helib/Ctxt.h>
#include <helib/ArgMap.h>

static long mValues[][14] = { 
//{ p, phi(m),  m,   d, m1, m2, m3,   g1,    g2,   g3,ord1,ord2,ord3, c_m}
  { 2,  512,    771, 16,771,  0,  0,     5,    0,    0,-32,  0,  0, 100}, // m=(3)*{257} :-( m/phim(m)=1.5 C=77 D=2 E=4
  { 2, 4096,   4369, 16, 17, 257, 0,   258, 4115,    0, 16,-16,  0, 100}, // m=17*(257) :-( m/phim(m)=1.06 C=61 D=3 E=4
  { 2, 16384, 21845, 16, 17, 5, 257,  8996,17477,21591, 16,  4,-16,1600}, // m=5*17*(257) :-( m/phim(m)=1.33 C=65 D=4 E=4
  { 2, 23040, 28679, 24, 17, 7, 241, 15184, 4098,28204, 16,  6,-10,1500}, // m=7*17*(241) m/phim(m)=1.24    C=63  D=4 E=3
  { 2, 46080, 53261, 24, 17,13, 241, 43863,28680,15913, 16, 12,-10, 100}, // m=13*17*(241) m/phim(m)=1.15   C=69  D=4 E=3
  { 2, 64512, 65281, 48, 97,673,  0, 43073,22214,    0, 96,-14,  0, 100}  // m=97*(673) :-( m/phim(m)=1.01  C=169 D=3 E=4 
};

#ifdef DEBUG_PRINTOUT
extern SecKey* dbgKey;
extern EncryptedArray* dbgEa;
#define FLAG_PRINT_ZZX  1
#define FLAG_PRINT_POLY 2
#define FLAG_PRINT_VEC  4
extern void decryptAndPrint(ostream& s, const Ctxt& ctxt, const SecKey& sk,
			    const EncryptedArray& ea, long flags=0);
#endif

std::string uint8_to_hex_string(const uint8_t *v, const size_t s);
void printState(Vec<uint8_t>& st);
extern long AESKeyExpansion(unsigned char RoundKey[],
			    unsigned char Key[], int NN);
extern void Cipher(unsigned char out[16],
		   unsigned char in[16], unsigned char RoundKey[], int Nr);

int main(int argc, char **argv)
{
  ArgMap amap;

  std::string plaintext; // PIIs of IDH
  amap.arg("pii", plaintext, "PIIs of IDH");

  std::string outfileName = "output.bin";
  amap.arg("outfileName", outfileName, "Output file for encrypted PIIs");

  long idx = 0;
  amap.arg("sz", idx, "parameter-sets: toy=0 through huge=5");

  long c=3;
  amap.arg("c", c, "number of columns in the key-switching matrices");

  long N=0;
  amap.arg("N", N, "# of bits of the modulus chain");

  long B=23;
  amap.arg("B", B, "# of bits per level (only 64-bit machines)");

  bool boot=false;
  amap.arg("boot", boot, "includes bootstrapping");

  bool packed=true;
  amap.arg("packed", packed, "use packed bootstrapping");

  amap.parse(argc, argv);

  if (idx>5) idx = 5;

  Vec<long> mvec;
  vector<long> gens;
  vector<long> ords;

  long p = mValues[idx][0];
  //  long phim = mValues[idx][1];
  long m = mValues[idx][2];

  append(mvec, mValues[idx][4]);
  if (mValues[idx][5]>1) append(mvec, mValues[idx][5]);
  if (mValues[idx][6]>1) append(mvec, mValues[idx][6]);

  gens.push_back(mValues[idx][7]);
  if (mValues[idx][8]>1)   gens.push_back(mValues[idx][8]);
  if (mValues[idx][9]>1) gens.push_back(mValues[idx][9]);

  ords.push_back(mValues[idx][10]);
  if (abs(mValues[idx][11])>1) ords.push_back(mValues[idx][11]);
  if (abs(mValues[idx][12])>1) ords.push_back(mValues[idx][12]);

  cout << "*** Test_AES: c=" << c
       << ", N=" << N
       << ", boot=" << boot
       << ", packed=" << packed
       << ", m=" << m
       << " (=" << mvec << "), gens="<<gens<<", ords="<<ords
       << endl;

  setTimersOn();
  double tm = -GetTime();
  cout << "computing key-independent tables..." << std::flush;
  Context context(m, p, /*r=*/1, gens, ords);
#if (NTL_SP_NBITS>=50) // 64-bit machines
  //context.bitsPerLevel = B;
#endif
  context.zMStar.set_cM(mValues[idx][13]/100.0); // the ring constant
  buildModChain(context, N, c);

  if (boot) context.makeBootstrappable(mvec);
  tm += GetTime();
  cout << "done in "<<tm<<" seconds\n";

  //  context.zMStar.printout();
  {IndexSet allPrimes(0,context.numPrimes()-1);
   cout <<"  "<<context.numPrimes()<<" primes ("
       <<context.ctxtPrimes.card()<<" ctxt/"
       <<context.specialPrimes.card()<<" special), total bitsize="
	<<context.logOfProduct(allPrimes)
	<<", security level: "<<context.securityLevel() << endl;}

  long e = mValues[idx][3] /8; // extension degree
  cout << "  "<<context.zMStar.getNSlots()<<" slots ("
       << (context.zMStar.getNSlots()/16)<<" blocks) per ctxt";
  if (boot && packed)
    cout << ". x"<<e<<" ctxts";
  cout << endl;

  cout << "computing key-dependent tables..." << std::flush;
  tm = -GetTime();
  SecKey secretKey(context);
  PubKey& publicKey = secretKey;
  secretKey.GenSecKey(64);      // A Hamming-weight-64 secret key

  // Add key-switching matrices for the automorphisms that we need
  long ord = context.zMStar.OrderOf(0);
  for (long i = 1; i < 16; i++) { // rotation along 1st dim by size i*ord/16
    long exp = i*ord/16;
    long val = PowerMod(context.zMStar.ZmStarGen(0), exp, m); // val = g^exp

    // From s(X^val) to s(X)
    secretKey.GenKeySWmatrix(1, val);
    if (!context.zMStar.SameOrd(0))
      // also from s(X^{1/val}) to s(X)
      secretKey.GenKeySWmatrix(1, InvMod(val,m));
  }

  addFrbMatrices(secretKey);      // Also add Frobenius key-switching
  if (boot) { // more tables
    addSome1DMatrices(secretKey);   // compute more key-switching matrices
    secretKey.genRecryptData();
  }
  tm += GetTime();
  cout << "done in "<<tm<<" seconds\n";

#ifdef DEBUG_PRINTOUT
  dbgKey = &secretKey; // debugging key and ea

  ZZX aesPoly;         // X^8+X^4+X^3+X+1
  SetCoeff(aesPoly,8);  SetCoeff(aesPoly,4);
  SetCoeff(aesPoly,3);  SetCoeff(aesPoly,1);  SetCoeff(aesPoly,0);
  dbgEa = new EncryptedArray(context, aesPoly);
#endif

  cout << "computing AES tables..." << std::flush;
  tm = -GetTime();
  HomAES hAES(context); // compute AES-specific key-independent tables
  const EncryptedArrayDerived<PA_GF2>& ea2 = hAES.getEA();
  long blocksPerCtxt = ea2.size() / 16;

  long nBlocks;
  if (boot && packed)
    nBlocks = blocksPerCtxt * e;
  else
    nBlocks = blocksPerCtxt;

  Vec<uint8_t> ptxt(INIT_SIZE, nBlocks*16);

  // Transform input plaintext 
  for(int i =0; i<ptxt.length(); i++){
    if(i < plaintext.length())
      ptxt[i] = plaintext[i];
    else
      ptxt[i] = 0;
  }


  // Encrypt the ptxt under the HE key
  vector< Ctxt > doublyEncrypted;
  tm = -GetTime();

  Vec<ZZX> encoded;
  encode4AES(encoded, ptxt, ea2);      // encode as HE plaintext

  Ctxt tmpCtxt(publicKey);
  cout << doublyEncrypted.size() << endl;
  doublyEncrypted.resize(encoded.length(), tmpCtxt); // allocate space
  cout << doublyEncrypted.size() << endl;
  for (long i=0; i<(long) doublyEncrypted.size(); i++) // encrypt the encoded key
    publicKey.Encrypt(doublyEncrypted[i], encoded[i]);


  // Write result

  std::ofstream outfile (outfileName);
  doublyEncrypted.data()->write(outfile);
  outfile.close();
  cout << "Done. Wrote result to " << outfileName << endl;

  std::ofstream publicKeyFile ("pubKey.bin");
  writePubKeyBinary(publicKeyFile, publicKey);
  publicKeyFile.close();

  std::ofstream secretKeyFile ("secretKey.bin");
  writeSecKeyBinary(secretKeyFile, secretKey);
  secretKeyFile.close();

  return 0;
}

#include <iomanip> 
void printState(Vec<uint8_t>& st)
{
  cerr << "[";
  for (long i=0; i<st.length() && i<32; i++) {
    cerr << std::hex << std::setw(2) << (long) st[i] << " ";
  }
  if (st.length()>32) cerr << "...";
  cerr << std::dec << "]";
}