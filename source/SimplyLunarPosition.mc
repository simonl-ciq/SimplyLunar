/*
Copyright(c) 2010, Darin Koblick
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met :
*Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using Toybox.Math;
using Toybox.System as Sys;

using Phase as phase;

const MYPI = 3.14159265358979323846d;

(:glance)
    module Astro {

        /**
        * php translation of a cart2sph MATLAB function
        * [azimuth,elevation,r] = cart2sph(x,y,z)
        * 
        * @author Dubravko Loborec
        * @copyright  2017, Dubravko Loborec
        * @link http://www.dubravkoloborec.com
        * 
        * @param mixed x
        * @param mixed y
        * @param mixed z
        * @return array[azimuth, elevation, r]
        */
        function cart2sph(x,y,z){
            var azimuth = Math.atan2(y,x);
            var elevation = Math.atan2(z, Math.sqrt(Math.pow(x,2) + Math.pow(y,2)));
            var r = Math.sqrt(Math.pow(x,2) + Math.pow(y,2) + Math.pow(z,2));
            return [azimuth, elevation, r];
        }

        /**
        * php translation of a sph2cart MATLAB function
        * [x,y,z] = sph2cart(azimuth,elevation,r)
        *  
        * @author Dubravko Loborec
        * @copyright  2017, Dubravko Loborec
        * @link http://www.dubravkoloborec.com
        * 
        * @param mixed azimuth
        * @param mixed elevation
        * @param mixed r
        * @return array[x, y, z]
        */
        function sph2cart(azimuth, elevation, r){
            var x = r * Math.cos(elevation) * Math.cos(azimuth);
            var y = r * Math.cos(elevation) * Math.sin(azimuth);
            var z = r * Math.sin(elevation);
            return [x, y, z];
        }

        /**
        * php translation of a dot MATLAB function
        *      
        * @author Dubravko Loborec
        * @copyright  2017, Dubravko Loborec
        * @link http://www.dubravkoloborec.com
        * 
        * @param array a
        * @param array b
        * @return float
        */
        function dot(a, b){
            var sum=0;
            for (var x = 0; x < a.size(); x++) {
                sum=sum+a[x]*b[x];
            } 
            return sum;   
        }

        /**
        * Returns integer value of the value
        * 
        * @param mixed value
        * @return int
        */
/*
        function bigintval(value){
			return Math.floor(value).toNumber();
        }
*/
        /**
        * helper function
        * 
        * @param mixed a
        * @param mixed b
        * @param mixed c
        */
        function sol_calc(a, b, c){
            b[0]=b[0].toDouble()*c.toDouble();  
            b[1]=b[1].toDouble()*c.toDouble();  
            b[2]=b[2].toDouble()*c.toDouble();
            var x = a[0][0]*b[0]+a[1][0]*b[1]+a[2][0]*b[2];
            var y = a[0][1]*b[0]+a[1][1]*b[1]+a[2][1]*b[2];
            var z = a[0][2]*b[0]+a[1][2]*b[1]+a[2][2]*b[2];

            return [x, y, z]; 
        }

        /**
        * Calculates juliandate from UTC date
        * php translation by Dubravko Loborec
        * 
        * @author Darin Koblick
        * @copyright  2010, Darin Koblick
        * @link https://www.mathworks.com/matlabcentral/fileexchange/22992-lunar-azimuth-and-altitude-estimation-algorithm
        * @license https://www.mathworks.com/matlabcentral/fileexchange/22992-lunar-azimuth-and-altitude-estimation-algorithm
        * 
        * @param integer y //UTC year
        * @param integer m //UTC month
        * @param integer d //UTC year
        * @param integer h //UTC hour
        * @param integer mn //UTC minute
        * @param integer s //UTC second
        * @return float
        */
//        function juliandate(y, m, d, h=0, mn=0, s=0){
        function juliandaynumber(y, m, d, h, mn, s) {
            if ( m <= 2 ){ // january & february
                y  = y - 1.0;
                m = m + 12.0;
            }
            return Math.floor(365.25d*(y + 4716.0d)).toNumber() + Math.floor(30.6001d*(m + 1.0d)).toNumber() + 2.0d - 
            Math.floor(y/100.0d).toNumber() + Math.floor(Math.floor(y/100.0d).toNumber()/4.0d).toNumber() + d - 1524.5d + 
            (h + mn/60.0d + s/3600.0d)/24.0d;
        }            

// % operator returns remainder a/b  (for php fmod())
		function myMod(a, b) {
			var d = (b < 0) ? -b : b;
			var m = (a - ((a / d).toLong() * d));
			var r = ((m < 0) ? (d + m) : m);

			return ((b < 0) ? (r + b) : r);
		}

		function array_sum(array) {
			var sum = 0;
			for (var i = 0; i < array.size(); i++) {
				sum += array[i];
			}
			return sum;
		}
        
        /**
        * Predict the Lunar Azimuth and Altitude within +/- .2 deg of any lat and lon for a given UTC
        * This algorithm will accept a Latitude, Longitude and Altitude location as well as a specific universal coordinated time. It will use this information and calculate the position of the moon in a local coordinate frame (az and alt aka az and el).
        * 
        * php translation and modification by Dubravko Loborec
        * 
        * @author Darin Koblick
        * @copyright  2010, Darin Koblick
        * @link https://www.mathworks.com/matlabcentral/fileexchange/22992-lunar-azimuth-and-altitude-estimation-algorithm
        * @license https://www.mathworks.com/matlabcentral/fileexchange/22992-lunar-azimuth-and-altitude-estimation-algorithm
        * 
        * @param integer y //UTC year
        * @param integer m //UTC month
        * @param integer d //UTC year
        * @param integer h //UTC hour
        * @param integer mn //UTC minute
        * @param integer s //UTC second
        * @param float Lat //Site Latitude in degrees -90:90 -> S(-) N(+) 
        * @param float Lon //Site Longitude in degrees -180:180 W(-) E(+) 
        * @param float Alt //Site Altitude in km
        * @return array[Az, h]
        */
//        function LunarAzEl(y, m, d, h=0, mn=0, s=0, Lat, Lon, Alt){
        function LunarAzEl(y, m, d, h, mn, s, Lat, Lon, Alt){

            while (Lon > 180){
                Lon = Lon - 360;
            }
            while (Lon < -180){
                Lon = Lon + 360; 
            }
            while (Lat > 90){
                Lat = Lat - 360;
            }
            while (Lat < -90){
                Lat = Lat + 360; 
            }

            var EarthRadEq = 6378.1370d;

            var jd = juliandaynumber(y, m, d, h, mn, s);
            var sd = jd - 2451543.5d;

            var N = 125.1228d-0.0529538083d*sd; //    (Long asc. node deg)
            var i = 5.1454d; //                      (Inclination deg)
            var w = 318.0634d + 0.1643573223d*sd; //  (Arg. of perigee deg)
            var a =  60.2666d;//                     (Mean distance (Earth's Equitorial Radii)
            var e = 0.054900d;//                     (Eccentricity)
            var M = myMod(115.3654d+13.0649929509d*sd,360);//    (Mean anomaly deg)

            var LMoon =  myMod(N + w + M,360);                 //(Moon's mean longitude deg)
            var FMoon =  myMod(LMoon - N,360);                 //(Moon's argument of latitude)

            //Keplerian Elements of the Sun
            var wSun = myMod(282.9404d + 4.70935E-5*sd,360);    // (longitude of perihelion)
            var MSun = myMod(356.0470d + 0.9856002585d*sd,360);  // (Sun mean anomaly)
            var LSun = myMod(wSun + MSun,360);                 // (Sun's mean longitude)

            var DMoon =  LMoon - LSun;                     // (Moon's mean elongation)  

            //Calculate Lunar perturbations in Longitude
            var LunarPLon = [ 
                -1.274*Math.sin((M - 2 *DMoon)*MYPI /180.0d), 
                .658*Math.sin(2 *DMoon*MYPI /180.0d), 
                -0.186*Math.sin(MSun*MYPI /180.0d), 
                -0.059*Math.sin((2 *M-2 *DMoon)*MYPI /180.0d), 
                -0.057*Math.sin((M-2*DMoon + MSun) *MYPI /180.0d), 
                .053 *Math.sin((M+2*DMoon) *MYPI /180.0d), 
                .046 *Math.sin((2*DMoon-MSun) *MYPI /180.0d), 
                .041 *Math.sin((M-MSun) *MYPI /180.0d), 
                -0.035 *Math.sin(DMoon *MYPI /180.0d),            
                -0.031 *Math.sin((M+MSun) *MYPI /180.0d), 
                -0.015 *Math.sin((2 *FMoon-2 *DMoon) *MYPI /180.0d), 
                .011 *Math.sin((M-4 *DMoon) *MYPI /180.0d)
            ];

            //Calculate Lunar perturbations in Latitude 
            var LunarPLat = [ 
                -0.173 *Math.sin((FMoon-2 * DMoon) *MYPI /180.0d), 
                -0.055 *Math.sin((M-FMoon-2 * DMoon) *MYPI /180.0d), 
                -0.046 *Math.sin((M+FMoon-2 * DMoon) *MYPI /180.0d), 
                +0.033 *Math.sin((FMoon+2 * DMoon) *MYPI /180.0d), 
                +0.017 *Math.sin((2 * M + FMoon) *MYPI /180.0d)
            ];

            //Calculate perturbations in Distance
            var LunarPDist = [                                   
                -0.58*Math.cos((M-2 * DMoon) * MYPI /180.0d),
                -0.46 * Math.cos(2 * DMoon * MYPI /180.0d)
            ];

            // Compute EE, the eccentric anomaly

            //E0 is the eccentric anomaly approximation estimate 
            //(this will initially have a relativly high error)
            var E0 = M+(180 / MYPI ) * e * Math.sin(M * MYPI /180.0d) * (1+e * Math.cos(M * MYPI /180.0d));

            //Compute E1 and set it to E0 until the E1 == E0
            var E1 = E0-(E0-(180 / MYPI ) * e * Math.sin(E0 * MYPI /180.0d)-M)/(1-e*Math.cos(E0 * MYPI /180.0d));

            while (E1-E0 > .000005){
                E0 = E1;
                E1 = E0-(E0-(180 / MYPI ) * e * Math.sin(E0 * MYPI /180.0d)-M)/(1-e*Math.cos(E0 * MYPI /180.0d));    
            }

            var EE = E1;

            //Compute rectangular coordinates (x,y) in the plane of the lunar orbit
            var x = a * (Math.cos(EE * MYPI /180.0d)-e);
            var sy = a * Math.sqrt(1-e * e) * Math.sin(EE * MYPI /180.0d);

            //convert this to distance and true anomaly
            var r = Math.sqrt(x * x + sy * sy);
            var v = Math.atan2(sy * MYPI /180.0d, x * MYPI /180.0d) * 180 / MYPI;

            //Compute moon's position in ecliptic coordinates
            var xeclip = r * (Math.cos(N * MYPI /180.0d) * Math.cos((v+w) * MYPI /180.0d) - Math.sin(N * MYPI /180.0d) * Math.sin((v+w) * MYPI /180.0d) * Math.cos(i * MYPI /180.0d));
            var yeclip = r * (Math.sin(N * MYPI /180.0d) * Math.cos((v+w) * MYPI /180.0d) + Math.cos(N * MYPI /180.0d) * Math.sin((v+w) * MYPI /180.0d) * Math.cos(i * MYPI /180.0d));
            var zeclip = r * Math.sin((v+w) * MYPI /180.0d) * Math.sin(i * MYPI /180.0d);

            //Add the calculated lunar perturbation terms to increase fmodel fidelity
            var list = cart2sph(xeclip,yeclip,zeclip);
            var eLon = list[0], eLat = list[1], eDist = list[2];

            list = sph2cart(eLon + array_sum(LunarPLon) * MYPI /180.0d, 
                eLat + array_sum(LunarPLat) * MYPI /180.0d, 
                eDist + array_sum(LunarPDist));
            xeclip = list[0]; yeclip = list[1]; zeclip = list[2];

            //  clear eLon eLat eDist;
            eLon=0;
            eLat=0;
            eDist=0;

            //convert the latitude and longitude to right ascension RA and declination
            //delta
            var T = (jd-2451545.0d)/36525.0d;

            //Generate a rotation matrix for ecliptic to equitorial
            //RotM=rotm_coo('e',jd);
            //See rotm_coo.M for obl and rotational matrix transformation
            var Obl = 23.439291d - 0.0130042d * T - 0.00000016d * T * T + 0.000000504d * T * T * T;
            Obl = Obl * MYPI /180.0d;

            //RotM = [1 0 0; 0 Math.cos(Obl) sin(Obl); 0 -sin(Obl) Math.cos(Obl)]'; /////////////////////

            var RotM= [[1, 0, 0],[ 0, Math.cos(Obl), Math.sin(Obl)],[0, -Math.sin(Obl), Math.cos(Obl)]];

            //Apply the rotational matrix to the ecliptic rectangular coordinates
            //Also, convert units to km instead of earth equatorial radii   
            var sol=sol_calc(RotM, [xeclip, yeclip, zeclip], EarthRadEq);
            //Find the equatorial rectangular coordinates of the location specified
            list = sph2cart(Lon * MYPI /180.0d,Lat * MYPI /180.0d,Alt+EarthRadEq);
            var xel = list[0], yel= list[1], zel = list[2];

            //Find the equatorial rectangular coordinates of the location @ sea level
            list = sph2cart(Lon * MYPI /180.0d,Lat * MYPI /180.0d,EarthRadEq);
            var xsl = list[0], ysl= list[1], zsl = list[2];

            //Find the Angle Between sea level coordinate vector and the moon vector
            /*theta1 = 180 - acosd(dot([xsl ysl zsl], [sol(1)-xsl sol(2)-ysl sol(3)-zsl]) 
            /(Math.sqrt(xsl.^2 + ysl.^2 + zsl.^2) 
            * Math.sqrt((sol(1)-xsl).^2 + (sol(2)-ysl).^2 + (sol(3)-zsl).^2)));*/

            var theta1 = 180.0d - Math.toDegrees(Math.acos(
                dot( [xsl, ysl, zsl], [sol[0]-xsl, sol[1]-ysl, sol[2]-zsl] )
                /(Math.sqrt(Math.pow(xsl,2) + Math.pow(ysl,2) + Math.pow(zsl,2)) 
                    * Math.sqrt( Math.pow(sol[0]-xsl,2) + Math.pow(sol[1]-ysl,2) + Math.pow(sol[2]-zsl,2)))

            ));         

            //Find the Angle Between the same coordinates but at the specified elevation
            /*theta2 = 180 - acosd(dot([xel yel zel],[sol(1)-xel sol(2)-yel sol(3)-zel]) ...
            ./(Math.sqrt(xel.^2 + yel.^2 + zel.^2) ...
            * Math.sqrt((sol(1)-xel).^2 + (sol(2)-yel).^2 + (sol(3)-zel).^2)));*/

            var theta2 = 180.0d - Math.toDegrees(Math.acos(
                dot( [xel, yel, zel], [sol[0]-xel, sol[1]-yel, sol[2]-zel] )
                /(Math.sqrt(Math.pow(xel,2) + Math.pow(yel,2) + Math.pow(zel,2)) 
                    *Math.sqrt( Math.pow(sol[0]-xel,2) + Math.pow(sol[1]-yel,2) + Math.pow(sol[2]-zel,2)))
            ));     

            //Find the Difference Between the two angles (+|-) is important
            var thetaDiff = theta2 - theta1;

            // equatorial to horizon coordinate tRAnsformation
            list = cart2sph(sol[0],sol[1],sol[2]);
            var RA = list[0], delta = list[1];
            delta = delta * (180 / MYPI );
            RA = RA * (180 / MYPI );

            //Following the RA DEC to Az Alt conversion sequence explained here:
            //http://www.stargazing.net/kepler/altaz.html

            //Find the J2000 value
            var J2000 = jd - 2451545.0d;

            //[Y,M,D,H,MN,S]
            //hourvec = datevec(UTC, 'yyyy/mm/dd HH:MM:SS');
            //UTH = hourvec(4) + hourvec(5)/60 + hourvec(6)/3600;

            var UTH = h + mn/60.0d + s/3600.0d;

//            var_dump(mn/60);
            //Calculate local siderial time
            var LST = myMod(100.46d+0.985647d * J2000+Lon+15*UTH,360);

            //Replace RA with hour angle HA
            var HA = LST-RA;

            //Find the h and Az at the current LST
            h = Math.asin(Math.sin(delta * MYPI /180.0d) * Math.sin(Lat * MYPI /180.0d) + Math.cos(delta * MYPI /180.0d) * Math.cos(Lat * MYPI /180.0d) * Math.cos(HA * MYPI /180.0d)) * (180 / MYPI );
            var Az = Math.acos((Math.sin(delta * MYPI /180.0d) - Math.sin(h * MYPI /180.0d) * Math.sin(Lat * MYPI /180.0d))/(Math.cos(h * MYPI /180.0d) * Math.cos(Lat * MYPI /180.0d))) * (180 / MYPI );

            //Add in the angle offset due to the specified site elevation
            h = h + thetaDiff;

            if (Math.sin(HA * MYPI /180.0d) >= 0){
                Az = 360.0d-Az; 
            }

            //Apply PaRAlax Correction if we are still on earth
            if (Alt < 100){
                var horParal = 8.794d/(r*6379.14d/149.59787e6);
                var p = Math.asin(Math.cos(h * MYPI /180.0d)*Math.sin((horParal/3600.0d) * MYPI /180.0d)) * (180 / MYPI );
                h = h-p;
            }
            return [Az, h];
        }

      

    }

