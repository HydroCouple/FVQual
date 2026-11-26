// SPDX-License-Identifier: MIT
// Small local reverse tape. Its lifetime is one solver operator/cell.
#pragma once
#include <cmath>
#include <vector>
#include <span>
#include <stdexcept>
namespace FVQual::Differentiation {
class Tape;
struct Scalar {
    double value=0;Tape *tape=nullptr;int index=-1;
    Scalar()=default;Scalar(double v):value(v){}
    explicit operator double() const {return value;}
};
struct Node {int a=-1,b=-1;double da=0,db=0;};
class Tape {
public:
    std::vector<Node> nodes;std::vector<int> inputs;
    Scalar input(double value) {Scalar x(value);x.tape=this;x.index=static_cast<int>(nodes.size());nodes.push_back({});inputs.push_back(x.index);return x;}
    static Scalar op(Scalar a,Scalar b,double value,double da,double db) {
        if(a.tape && b.tape && a.tape!=b.tape)throw std::invalid_argument("mixed local differentiation tapes");
        Scalar result(value);result.tape=a.tape?a.tape:b.tape;
        if(result.tape) {
            if(!std::isfinite(value)||!std::isfinite(da)||!std::isfinite(db))throw std::invalid_argument("nonfinite local derivative");
            result.index=static_cast<int>(result.tape->nodes.size());
            result.tape->nodes.push_back({a.index,b.index,da,db});
        }
        return result;
    }
    std::vector<double> jvp(std::span<const double> direction,std::span<const Scalar> outputs) const {
        if(direction.size()!=inputs.size())throw std::invalid_argument("local tangent dimensions differ");
        std::vector<double> d(nodes.size(),0),result;
        for(size_t k=0;k<inputs.size();++k)d[inputs[k]]=direction[k];
        for(size_t k=0;k<nodes.size();++k) {const auto &n=nodes[k];if(n.a>=0)d[k]+=n.da*d[n.a];if(n.b>=0)d[k]+=n.db*d[n.b];}
        for(const auto &x:outputs)result.push_back(x.index<0?0:d[x.index]);return result;
    }
    std::vector<double> vjp(std::span<const Scalar> outputs,std::span<const double> cotangents) const {
        if(outputs.size()!=cotangents.size())throw std::invalid_argument("local cotangent dimensions differ");
        std::vector<double> d(nodes.size(),0),result;
        for(size_t k=0;k<outputs.size();++k)if(outputs[k].index>=0)d[outputs[k].index]+=cotangents[k];
        for(size_t k=nodes.size();k>0;--k) {const auto &n=nodes[k-1];if(n.a>=0)d[n.a]+=n.da*d[k-1];if(n.b>=0)d[n.b]+=n.db*d[k-1];}
        for(int index:inputs)result.push_back(d[index]);return result;
    }
};
inline Scalar operator+(Scalar a,Scalar b){return Tape::op(a,b,a.value+b.value,1,1);}
inline Scalar operator-(Scalar a,Scalar b){return Tape::op(a,b,a.value-b.value,1,-1);}
inline Scalar operator*(Scalar a,Scalar b){return Tape::op(a,b,a.value*b.value,b.value,a.value);}
inline Scalar operator/(Scalar a,Scalar b){return Tape::op(a,b,a.value/b.value,1/b.value,-a.value/(b.value*b.value));}
inline Scalar operator-(Scalar a){return Scalar(0)-a;}
inline Scalar exp(Scalar a){const double y=std::exp(a.value);return Tape::op(a,{},y,y,0);}
inline Scalar pow(Scalar a,double b){return Tape::op(a,{},std::pow(a.value,b),b==0?0:b*std::pow(a.value,b-1),0);}
inline Scalar sin(Scalar a){return Tape::op(a,{},std::sin(a.value),std::cos(a.value),0);}
inline Scalar cos(Scalar a){return Tape::op(a,{},std::cos(a.value),-std::sin(a.value),0);}
inline Scalar log(Scalar a){return Tape::op(a,{},std::log(a.value),1/a.value,0);}
inline Scalar abs(Scalar a){return a.value<0?-a:a;}
inline Scalar pow(Scalar a,Scalar b){
    if(b.index<0)return pow(a,b.value);
    if(a.value<=0)throw std::invalid_argument("variable exponent requires a positive base for differentiation");
    const double y=std::pow(a.value,b.value);
    return Tape::op(a,b,y,b.value*y/a.value,y*std::log(a.value));
}
}
