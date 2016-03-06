//
//  main.cpp
//  calculateDewT
//
//  Created by Morgan on 15/11/8.
//  Copyright © 2015年 Morgan. All rights reserved.
//

//开尔文温度 -273.15℃=0K
//当前压力：currentPressure
//水的饱和温度
#include <iostream>
#include <math.h>
using namespace std;

//摄氏度转开尔文温度
double conversionT_C2K(double temperature_C);

//开尔文温度转摄氏度
double conversionT_K2C(double temperature_K);

//压力转换 mm汞柱转千帕
double conversion_P_lgmmHg2kPa(double lgp);

//计算水的饱和蒸汽压力 
double lgsvp(double saturationTemperatureH2O);

//计算水的饱和温度
double saturationTemperatureH2O(double saturationPressureH2O_kPa);

//千卡换算千焦
double conversionQ_kCal2kJ(double kCal);

//水的焓值（kJ/kg）
double _H2O_enthalpy(double temperatureH2O_C);

//水蒸气的焓值
double _H2OVapor_enthalpy(double temperatureH2O_C);

//过热水蒸汽的焓值
double _H2OHeat_enthalpy(double saturationTemperatureH2O_C, double saturationTemperatureLiBr_C);

//溴化锂溶液露点温度
double dewTLiBr(double saturationTemperatureH2O_C, double concentrationOfLiBrSolution);

//溴化锂水溶液焓值
double enthalpyLiBrSolution(double solutionTemperatureLiBr_C,double concentrationOfLiBrSolution);

//溴化锂溶液浓度
double _concentration_LiBrSolution(double solutionTemperatureLiBr_C,double pressure_kPa);

int main(int argc, const char * argv[]) {
    // insert code here...
    
    double currentPressure;
    double saturationTemperatureH2O_K;
    
    double saturationTemperatureH2O_C;
    double concentrationOfLiBrSolution;
    double dewTOfLiBr;
    
    double solutionTemperatureLiBr_C;
    double h_LibrSolution;
    double h_H2O;
    
//--求水的饱和蒸汽压力
//    cout << "请输入当前温度（摄氏度）\n";
//    cin >> saturationTemperatureH2O_C;
//    
//    double lgp;
//    saturationTemperatureH2O_K =  conversionT_C2K(saturationTemperatureH2O_C);
//    lgp = lgsvp(saturationTemperatureH2O_K);
//    currentPressure = conversion_P_lgmmHg2kPa(lgp);
//    
//    cout << "当前温度下水的饱和蒸汽压为：\n"
//    << "lgp =" <<lgp <<endl
//    <<"p ="<< currentPressure <<" kPa"<< endl;

//--求水的饱和蒸汽温度
//    cout << "请输入当前压力（kPA）\n";
//    cin >> currentPressure;
//    cout << "饱和温度T =" << saturationTemperatureH2O(currentPressure) << "˚C\n";
    
//--求溴化锂溶液的露点温度
    cout << "\n输入压力为p时水的饱和温度（摄氏度，0℃ < t < 100℃）\n";
    cin >> saturationTemperatureH2O_C;
    
    cout << "输入溴化锂水溶液中含有溴化锂的千克数（45% < x < 65%）\n";
    cin >> concentrationOfLiBrSolution;
    dewTOfLiBr = dewTLiBr(saturationTemperatureH2O_C,concentrationOfLiBrSolution);
    
    cout << "溴化锂溶液的露点温度为：\n";
    cout << "t =" << dewTOfLiBr<<"℃"<<endl;

//--求水和水蒸气的焓值
    cout << "此温度下水的焓值为: \n";
    cout << "h1 =" << _H2O_enthalpy(saturationTemperatureH2O_C) << "kJ/kg" <<endl;
    cout << "水蒸气的焓值为： \n";
    cout << "h2 =" << _H2OVapor_enthalpy(saturationTemperatureH2O_C) << "kJ/kg" << endl;
    cout << "过热水蒸气的焓值为： \n";
    cout << "h = " << _H2OHeat_enthalpy(128.744, 75);
    
//--求溴化锂溶液的浓度
//    cout << "\n输入溴化锂溶液的温度(˚C)\n";
//    cin >> solutionTemperatureLiBr_C;
//    cout << "输入溴化锂溶液的压强(kPa)\n";
//    cin >> currentPressure;
//    double x = _concentration_LiBrSolution(solutionTemperatureLiBr_C, currentPressure);
//    cout  << "x =" << x << "%" <<endl;
//    
//--求溴化锂溶液的焓值
//    cout << "输入溴化锂溶液的温度(˚C)\n";
//    cin >> solutionTemperatureLiBr_C;
//    
//    cout << "输入溴化锂水溶液中含有溴化锂的千克数（30% < x < 75%）\n";
//    cin >> concentrationOfLiBrSolution;
//    
//    h_LibrSolution = enthalpyLiBrSolution(solutionTemperatureLiBr_C, concentrationOfLiBrSolution);
//    
//    cout << "溴化锂溶液的焓值为：\n";
//    cout << "h =" << h_LibrSolution << " kj/kg" <<endl;
    

    
    return 0;
}

//END MAIN
//----------------------------------------------------------------------------------

/*
 #01
 摄氏度转开尔文温度
 */
double conversionT_C2K(double temperature_C){
    return temperature_C + 273.15;
}

/*
 #02
 开尔文温度转摄氏度
 */
double conversionT_K2C(double temperature_K){
    return temperature_K - 273.15;
}

/*
 #03
 mmHg与kPa换算
 */
double conversion_P_lgmmHg2kPa(double lgp){
    return  pow(10, lgp)*0.1333224;
}

/*
 #04
 计算水的饱和蒸汽压力  公式来源 硕士论文《溴化锂热泵设计》- 杨玉华  （1）
  Calculate saturated vapor pressure(svp) of H2O.
  lgsvp，lg为以10为底的对数
*/
double lgsvp(double saturationTemperatureH2O_K){
    double lgsvp;
    double t = saturationTemperatureH2O_K;
    double lgT = log10(t);
    lgsvp = 31.46559 - 8.2*lgT - 3142.305/t + 0.0024804*t;
    return lgsvp;
}

/*
 计算水的饱和蒸汽压力  公式来源 本科生毕业论文《第一类溴化锂吸收式热泵设计.doc》
 Calculate saturated vapor pressure(svp) of H2O (˚C).
 返回kPa
 此公式来误差太大 t=27 时 p=2906.39kpa.
 */
//double svp(double saturationTemperatureH2O_C){
//    double p;
//    double t = saturationTemperatureH2O_C;
//    double a = 9.4865;
//    double b = 389.27/(42.6676-t-273);
//    p=exp(a+b);
//    return p;
//}


/*计算水的饱和温度
 Calculate saturated vapor temperature(svp) of H2O.
 适用于0.1Mpa- 0.35Mpa 相对误差≤0.3%
 */
//double saturationTemperatureH2O(double saturationPressureH2O_kPa){
//    double t;
//    double p = saturationPressureH2O_kPa/1000;
//    
//    t = 189.2000 * sqrt(sqrt(p)) + -6.5;
//
//    return t;
//}

/*
 #05
 计算水的饱和温度 公式来源《蒸发过程饱和水蒸气对应参数的计算方程》-高俊
 Calculate saturated vapor temperature(svp) of H2O.
 适用于2kpa- 1500kpa 相对误差≤0.3%
 按此公式与上述压力计算公式反算 每个温度值与标准值几乎都相差0.4，遂在计算末加之
 */
double saturationTemperatureH2O(double saturationPressureH2O_kPa){
    double t = 0.0000000;
    double p = saturationPressureH2O_kPa;
    double c[] = {6.004553,15.89247,-0.1723261,0.3940772,-4.631586e-2,3.01707e-3};
    for (int i = 0; i < 6; i++) {
        t += c[i]*pow(log(p), i);
    }
    
    return t+0.4;
}

/*
 #06
 换算：1千卡 = 4.1858千焦
 */
double conversionQ_kCal2kJ(double kCal){
    return kCal * 4.1858;
}


/*
 #07
 计算饱和水的气化潜热(latent heat of vaporization)
 r = 597.34 - 0.555t1 - 0.2389 * 10^y
 y = 5.1463 - 1540/(t1+273.16)
 结果转化为kJ
 t1 为压力p时饱和水蒸气的温度，可利用#05计算
 */
double _H2O_latentHeatOfVaporization(double saturationTemperatureH2O_C){
    
    double t1 = saturationTemperatureH2O_C;
    double y = 5.1463 - 1540/(t1+273.16);
    double r = (597.34 - 0.555 * t1 - 0.2389 * pow(10, y))*4.1858;
//    cout << "\t温度为" << t1 << "时水的气化潜热r = " << r << endl;
    return r;
}

/*
 #08
 计算水的焓值: h1 = t1 + 100 (kcal/kg)
 t1为水的温度（˚C）
 结果用公式6换算成千焦
*/
double _H2O_enthalpy(double temperatureH2O_C){
    return conversionQ_kCal2kJ(temperatureH2O_C + 100.00);
}

/*
 #09
 计算水蒸气的焓值：h2 = h1 + r
 r为饱和水的气化潜热
*/
double _H2OVapor_enthalpy(double temperatureH2O_C){
    return _H2O_enthalpy(temperatureH2O_C) + _H2O_latentHeatOfVaporization(temperatureH2O_C);
}

/*
 #0
过热水蒸气的焓值
 h = h2 + Cp(t-t1)
 h2 = h1 + r
 h1 = t1 + 100
 
 t1 为压力p时饱和水蒸气的温度，可利用#05计算 ˚C
 t 为过热水蒸气的温度即压力p时溶液的平衡温度，可利用#08计算 ˚C
 h1 为温度为t1时饱和水的焓 kcal/kg
 h2 为温度为kcal/kg
 Cp 为过热水蒸气t1-t的定呀平均比热。计算普通单级循环蒸汽焓值时,Cp ≈ 0.46kcal/kg
 */

double _H2OHeat_enthalpy(double saturationTemperatureH2O_C, double saturationTemperatureLiBr_C){
    double Cp = 0.46;
    double t = saturationTemperatureLiBr_C;
    double t1 = saturationTemperatureH2O_C;
    double h2 = _H2OVapor_enthalpy(saturationTemperatureH2O_C);
    double h = h2 + conversionQ_kCal2kJ(Cp * (t - t1));
    return h;
}



/*
 #
 计算溴化锂的露点温度(Calculate dew temperature of LiBr)

    saturationTemperatureH2O_C:水的饱和温度（摄氏度）
    concentrationOfLiBrSolution :100kg 溴化锂水溶液中含有溴化锂的千克数
 此方程试用范围：
    0℃ < saturationTemperatureH2O_C < 100℃
    45% < concentrationOfLiBrSolution < 65%
    
    x = 61.7, t = 75 ; dewTLiBr = 17.64793616
 
*/
double dewTLiBr(double saturationTemperatureH2O_C, double concentrationOfLiBrSolution){
    
    double a[]={0.770033,1.45455e-2,-2.63906e-4,2.27609e-6};
    double b[]={140.877,-8.55749,0.16709,-8.82641e-4};
    
//      double a[]={0.77,1.455,-2.6401,2.277};
//      double b[]={140.876,-855.745,1670.89,882.636};
    
    double t= saturationTemperatureH2O_C;
    double x = concentrationOfLiBrSolution;
    
    double sum1 =0.000000 ,sum2=0.000000;
    
    //循环语句用于求和
    for (int i =0; i<4; i++) {
        sum1+=t*a[i]*pow(x, i);
        sum2+=b[i]*pow(x, i);
    }

    double dewTLiBr = sum1 + sum2;
    
    return dewTLiBr;
}


/*计算溴化锂水溶液的焓值（Calculate the enthalpy of LiBr solution）单位：kcal/kg;
 
 solutionTemperatureLiBr_C: 溴化锂溶液的温度，˚C
 concentrationOfLiBrSolution :100kg 溴化锂水溶液中含有溴化锂的千克数
 
 45% < concentrationOfLiBrSolution < 65%
 
 1千卡(kcal)=1大卡=4.184千焦(kJ)
 */
//double enthalpyLiBrSolution(double solutionTemperatureLiBr_C,double concentrationOfLiBrSolution){
//    double a[] = {-121.189,   16.7809,     -0.517766,  6.34755e-3,  -2.60914e-4};
//    double b[] = {0.671458,   1.01548e-2,  5.41941e-4, 6.82514e-6,  -2.80048e-8};
//    double c[] = {1.23744e-3, -7.74557e-5, 1.94305e-6, 6.52880e-11, 6.52880e-11};
//    double t = solutionTemperatureLiBr_C;
//    double x = concentrationOfLiBrSolution;
//    double sum1 = 0.000000, sum2 = 0.000000,sum3 = 0.000000;
//    double h;
//    
//    //h = ∑Anx^n + t∑Bnx^n + t^2*∑Cnx^n; (n = 0,1,2,3,4)
//    for (int i = 0; i < 5; i++) {
//        sum1 += a[i]*pow(x, i);
//        sum2 += t*(b[i]*pow(x, i));
//        sum3 += t*t*(c[i]*pow(x,i));
//        
//        cout << sum1 << " " << sum2 << " " <<sum3 <<endl;
//        
//    }
//        cout << sum1 << " " << sum2 << " " <<sum3 <<endl;
//    
//    cout << a[4]*pow(x, 4) <<endl;
//    
//    h = sum1 + sum2 + sum3;
//    
//    return h;
//}
/*
    0˚C ≤ solutionTemperatureLiBr_C  ≤ 200˚C  
    30% ≤ concentrationOfLiBrSolution ≤ 75%
 */
double enthalpyLiBrSolution(double solutionTemperatureLiBr_C,double concentrationOfLiBrSolution){
    double a[] = {3.22313e2,3.83413e2,-2.65438e3,2.87262e3};
    double b[] = {4.19928,-9.39005,1.60770e1,-1.36171e1};
    double c[] = {1.00479e-3,-1.41857e-3,-2.06186e-3,5.92438e-3};
    double sum1 = 0.00000000 , sum2 = 0.0000000 , sum3 = 0.0000000;
    double t = solutionTemperatureLiBr_C;
    double x = concentrationOfLiBrSolution/100;
    double h;
    /*
     sum1 = ∑ax^n,
     t = 45.83
     x = 58.998
     h = 280.218
     */
//    for (int i=0; i<4; i++) {
//        sum1 += a[i]*pow(x,i);
//        sum2 += b[i]*pow(x, i);
//        if (i<3) {
//            sum3 += c[i]*pow(x,i);
//        }
//    }
//    
//    h = sum1 + sum2*t  + (sum3 + pow(c[3],3)*x)*t*t;
    
    /*
     sum1 = ∑ax^n,
     t = 45.83
     x = 58.998
     h = 282.773
     */
    
    for (int i=0; i<4; i++) {
        sum1 += a[i]*pow(x,i);
        sum2 += b[i]*pow(x, i);
        sum3 += c[i]*pow(x, i);
    }
    
    h = sum1 + sum2*t + sum3*t*t;
    
    return  h;
    
}

/*
 已知溶液的温度t和压强P,确定溶液的浓度X%,平均相对误差0.76%，最大相对误差1.749%
 t----压强为P时，溴化锂溶液的饱和温度，˚C
 t1----压强为P时，对应水的饱和温度，˚C，即需要先根据P求出水的饱和温度
 10˚ ≤ t ≤ 130˚C
 2kpa ≤ P ≤ 1500kpa
 */
double _concentration_LiBrSolution(double solutionTemperatureLiBr_C,double pressure_kPa){
    double t = solutionTemperatureLiBr_C;
    double t1 = saturationTemperatureH2O(pressure_kPa);
    double a[] = {0.31057,-1.282e-2,-1.7312e-4,5.3303e-7};
    double b[] = {1.232e-2,3.846e-4,-7.1457e-8,-5.73e-9};
    double c[] = {-1.9166e-4,-3.334e-6,5.3123e-8,-3.6012e-10,1.0257e-12};
    double d[] = {1.6386e-6,-2.16e-8,1.505e-10,-4.678e-13};
    double sum1 = 0.000000,sum2= 0.000000,sum3 = c[4]*pow(t1, 4),sum4 = 0.000000;
    
    double x = 0.0000000;
    
    for (int i = 0; i < 4; i++) {
        sum1 += a[i] * pow(t1, i);
        sum2 += b[i] * pow(t1, i);
        sum3 += c[i] * pow(t1, i);
        sum4 += d[i] * pow(t1, i);
    }
    
    x = sum1 + sum2*t + sum3*t*t + sum4*t*t*t;
    
    return x * 100;
}

