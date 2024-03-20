
#include <ArduinoOTA.h>
#ifdef ESP32
  #include <WiFi.h>
  #include <AsyncTCP.h>
#else
  #include <ESP8266WiFi.h>
  #include <ESPAsyncTCP.h>
#endif
#include <ESPAsyncWebServer.h>
#include"matfun.h"
#include"geofun.h"

AsyncWebServer server(80);

const char* ssid = "srihari";
const char* password = "12345678";

const char* input_parametera = "inputa";
const char* input_parameterb = "inputb";
const char* input_parameterc = "inputc";
const char* input_parameterd = "inputd";
const char* input_parametere = "inpute";
const char* input_parameterf = "inputf";

const char index_html[] PROGMEM = R"rawliteral(
<!DOCTYPE HTML>
<html>
<head>
  <title>MAT ALG</title>
  <meta name="viewport" content="width=device-width, initial-scale=2">
  <style>
   html{font-family:Serif; text-align:left}
   h2{font-size:30px;color:#0000ff}
   input{background-color:Aqua,margin:4px;}

  </style>
 <script>
  function updateResults(results){
   document.getElementById("results").innerHTML = results;
   }
  </script>
  
</head>
<body>
  <h2>Triangle Properities of Mat_alg</h2>
  <p>Input points A, B, C</p>
  <form id="triangleData" onsubmit="submitForm(); return false;">
    Enter the values of Point A: <input type="number" value="1" name="inputa"><input type="number" value = "-1" name="inputb"><br>
    Enter the values of Point B: <input type="number" value="-4" name="inputc"><input type="number" value="6" name="inputd"><br>
    Enter the values of Point C: <input type="number" value="-3" name="inpute"><input type="number" value="-5" name="inputf"><br>
    <input type="submit" value="Submit">
  </form><br>
  <div id="results"></div>
<script>


   function submitForm(){
    var formData = new FormData(document.getElementById("triangleData"));
    fetch('/get',{method: "POST", body: formData})
     .then(response => response.text())
     .then(results => updateResults(results));
     }
  </script>
  
</body>
</html>
)rawliteral";


void notFound(AsyncWebServerRequest *request) {
  request->send(404, "text/plain", "Not found");
}

void setup() {
  Serial.begin(115200);
  WiFi.mode(WIFI_STA);
  WiFi.begin(ssid, password);
  if (WiFi.waitForConnectResult() != WL_CONNECTED) {
    Serial.println("Connecting...");
    return;
  }
  Serial.println();
  Serial.print("IP Address: ");
  Serial.println(WiFi.localIP());

  server.on("/", HTTP_GET, [](AsyncWebServerRequest *request){
    request->send_P(200, "text/html", index_html);
  });

  server.on("/get", HTTP_POST, [] (AsyncWebServerRequest *request) {
    double X1 = request->arg(input_parametera).toDouble();
    double Y1 = request->arg(input_parameterb).toDouble();
    double X2 = request->arg(input_parameterc).toDouble();
    double Y2 = request->arg(input_parameterd).toDouble();
    double X3 = request->arg(input_parametere).toDouble();
    double Y3 = request->arg(input_parameterf).toDouble();


    
    int m=2,n=1,i,j,k=0;
    
    //Initializing matrices
    
    double **A = createMat(m,n);
    
    A[0][0] = X1;
    A[1][0] = Y1;
    
    double **B = createMat(m,n);
    B[0][0] = X2;
    B[1][0] = Y2;
    
    double **C = createMat(m,n);
    C[0][0] = X3;
    C[1][0] = Y3;
    
    
    //Medians 
 double **D = Matsec(A, B, 2, 1.0);
 double **E = Matsec(B, C, 2, 1.0);
 double **F = Matsec(C, A, 2, 1.0);           

 
 double **side_ab = Matsub(A,D,m,n);
 double **side_bc= Matsub(B,E,m,n);
 double **side_ca = Matsub(C,F,m,n);
 
 //Normalizing Matrices
 
 double sideAB = Matnorm(side_ab,m);
 double sideBC = Matnorm(side_bc,m); 
 double sideCA = Matnorm(side_ca,m);
 
 //Orthocenter
    double **m1=Matsub(A,B,m,n);
    double **m2=Matsub(B,C,m,n);
    double **m3=Matsub(C,A,m,n);
    
    double **p= createMat(2,1);
      
    double **N=transposeMat(combineMat(m3,m1,m,n),m,m);
    double **m3_t=transposeMat(m3,m,n);
    double **m1_t=transposeMat(m1,m,n);
    double **p0=Matmul(m3_t,B,n,m,n);
    double **p1=Matmul(m1_t,C,n,m,n);
    p[0][0]=p0[0][0];
    p[1][0]=p1[0][0];
    

    double **l = createMat(m,m);
    double **O = createMat(2,1);


    O[1][0]=(-1)*(p[0][0]*N[1][0]-N[0][0]*p[1][0])/Matdet(N);
    O[0][0]=(p[0][0]*N[1][1]-N[0][1]*p[1][0])/Matdet(N);
    
     //Perpedicular bisector   
 
 double **_ab = Matsub(A,B,m,n);
 double **_bc = Matsub(B,C,m,n);
 double **_ca = Matsub(C,A,m,n);
 double **bisectorABM = normVec(_ab);
 double **bisectorBCM = normVec(_bc);
 double **bisectorCAM= normVec(_ca);
 
 //Angular Bisector
 
 int ma1 = 2, na1 = 1;
 double **Inc = createMat(ma1, na1);
 double **diff_AB = Matsub(B, A, 2, 1); 
 double distance_AB = Matnorm(diff_AB, 2);  
 double **diff_AC = Matsub(C, A, 2, 1);
 double distance_AC = Matnorm(diff_AC, 2);  
 double **diff_BC = Matsub(C, B, 2, 1);
 double distance_BC = Matnorm(diff_BC, 2);

 
 double a = distance_BC;
 double b = distance_AC;
 double c = distance_AB;
 
 double l1 = (a + c - b) / 2;
 double l2 = (a + b - c) / 2;
 double l3 = (c+ b - a) / 2;
 
 
 double **T1 = Matscale(A, 2, 1, a);
 double **T2= Matscale(B, 2, 1, b);
 double **T3 = Matscale(C, 2, 1, c);
 
 
 double **eigenvalues = Matadd(Matadd(T1, T2, 2, 1), T3, 2, 1);
 double eigendenominator = a + b + c;
 
 //Incenter 
 
 double **Da, **Ea, **Fa;
 Da = Matscale(Matadd(Matscale(C, 2, 1, l1), Matscale(B, 2, 1, l2), 2, 1), 2, 1, 1.0 / (l1 + l2));
 Ea = Matscale(Matadd(Matscale(A, 2, 1, l2), Matscale(C, 2, 1, l3), 2, 1), 2, 1, 1.0 / (l2 + l3));
 Fa = Matscale(Matadd(Matscale(B, 2, 1, l3), Matscale(A, 2, 1, l1), 2, 1), 2, 1, 1.0 / (l3 + l1));

 Inc[0][0] = eigenvalues[0][0] / eigendenominator;
 Inc[1][0] = eigenvalues[1][0] / eigendenominator;
 
    // Send the results to client
    
    String result = " Median AD: " + String(sideAB, 2) + "<br>"
                     "Median BE: " + String(sideBC, 2) + "<br>"
                     "Median CF: " + String(sideCA, 2) + "<br>"
    
                      "Orthocenter O: " + String(O[0][0],2) + ", " + String(O[1][0], 2) + "<br>"
    
                      "Bisector AB: " + String(bisectorABM[0][0], 2) +", "+ String(bisectorABM[1][0], 2) + "<br>"
                      "Bisector BC: " + String(bisectorBCM[0][0], 2) +", "+ String(bisectorBCM[0][0], 2)+ "<br>"
                      "Bisector CA: " + String(bisectorCAM[0][0], 2) +", "+ String(bisectorCAM[0][0], 2) + "<br>"
                      
   
                      "Bisector of A :" + String(D[0][0],2) + ", " + String(D[1][0],2) + "<br>"
                      "Bisector of B :" + String(E[0][0],2) + ", " + String(E[1][0],2) + "<br>"
                      "Bisector of C:" + String(F[0][0],2) + ", " + String(F[1][0],2) + "<br>"
                      
                      
                      "Distance of AB: " + c + "<br>"
                      "Distance of BC: " + a + "<br>"
                      "Distance of CA: " + b + "<br>"
                      
                      
                      "InCenter: "  + String(Inc[0][0],2) + ", " + String(Inc[1][0], 2) + "<br>";
    
    request->send(200, "text/html", result);
  });

  server.onNotFound(notFound);
  server.begin();
}

void loop() {
  // Nothing to do here
}

