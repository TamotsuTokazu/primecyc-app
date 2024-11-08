#include "binfhecontext.h"
#include "math/math-hal.h"

using Vector = lbcrypto::BigVector;
using Integer = Vector::Integer;

const Integer t("65");
const usint n = 100;
const usint p0 = 257;
const usint p1 = 433;
const usint pq = p0 * p1;
const usint Bks = 1 << 6;
const Integer Q("72057594422075861");


const Integer rootOfUnitypq("56503092529895247");
const Integer rootOfUnity0("17774520565003034");
const Integer rootOfUnity1("71367447949685206");


const Integer bigModulus("1393796574908163946345982392040522604871681");
const Integer bigRootOfUnity0("742774059212952446954038657883272672063730");
const Integer bigRootOfUnity1("742774059212952446954038657883272672063730");
const Integer bigRootOfUnitypq("1296198229327671723145458703434013327936548");

// const Integer t("65");
// const usint n = 1;
// const usint p0 = 5;
// const usint p1 = 7;
// const usint pq = p0 * p1;
// const usint Bks = 1 << 6;
// const Integer Q("72057594037958621");


// const Integer rootOfUnitypq("4380087458826051");
// const Integer rootOfUnity0("51233122604688005");
// const Integer rootOfUnity1("10890886932265458");


// const Integer bigModulus("680564733841876926926749214863536426241");
// const Integer bigRootOfUnity0("297199282715736452036733242311617415102");
// const Integer bigRootOfUnity1("297199282715736452036733242311617415102");
// const Integer bigRootOfUnitypq("584311738605810263590153731555488261050");

int main() {
    Vector a((p0 - 1) * (p1 - 1), Q);

    a[0] = 1;

    ChineseRemainderTransformArb<Vector>().SetCylotomicPolynomial(lbcrypto::GetCyclotomicPolynomial<Vector>(pq, Q), Q);

    Vector b = ChineseRemainderTransformArb<Vector>().ForwardTransform(a, rootOfUnitypq, bigModulus, bigRootOfUnitypq, pq);
    Vector c = ChineseRemainderTransformArb<Vector>().InverseTransform(b, rootOfUnitypq, bigModulus, bigRootOfUnitypq, pq);

    // std::cout << "a: " << a << std::endl;
    // std::cout << "b: " << b << std::endl;
    std::cout << "c: " << c << std::endl;

    return 0;
}