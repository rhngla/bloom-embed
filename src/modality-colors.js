// Shared count scale spans the observed RNA and ATAC endpoint values.
const umiValues=D.cells.flatMap(c=>[c[5],c[6]]),umiMin=Math.min(...umiValues),umiMax=Math.max(...umiValues);
const umiStops=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[230,205,40]];
let umiScale='linear';
function ease(v){return v*v*(3-2*v)}
function paletteColor(u){const v=Math.max(0,Math.min(1,u))*4,k=Math.min(3,Math.floor(v)),f=v-k;return `rgb(${umiStops[k].map((c,i)=>Math.round(c+(umiStops[k+1][i]-c)*f)).join(',')})`}
const umiPalette=Array.from({length:256},(_,i)=>paletteColor(i/255));
function countExtent(){return umiScale==='log'?[Math.log(umiMin),Math.log(umiMax)]:[umiMin,umiMax]}
function blendedCount(i,v=t){const c=D.cells[i],u=ease(v),a=umiScale==='log'?Math.log(c[5]):c[5],b=umiScale==='log'?Math.log(c[6]):c[6];return a+(b-a)*u}
function umiColor(i){const [lo,hi]=countExtent();return umiPalette[Math.round((blendedCount(i)-lo)/(hi-lo||1)*255)]}
function formatCount(v){return umiScale==='log'?v.toFixed(2):Math.round(v).toLocaleString()}
