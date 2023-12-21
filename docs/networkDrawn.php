<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<title>RecoPhyNC - Recognition of Phylogenetic Network Classes</title>
  <meta name="Description" content="RecoPhyNC - Recognition of Phylogenetic Network Classes">
  <link rel="stylesheet" type="text/css" media="screen" href="style.css">
  <link rel="shortcut icon" type="image/x-icon" href="favicon.ico">
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8"></head>

<body style="font-family:Calibri,Arial,sans-serif">

<script type="text/vnd.graphviz" id="cluster">

<?php
echo "digraph G{";

if(isset($_POST['DOTcode'])){
   $dotCode=$_POST['DOTcode'];


$text = trim($dotCode);
$textAr = explode("\n", $text);
$textAr = str_replace("\r","",$textAr);
//$textAr = array_filter($textAr, 'trim'); // remove any extra \r characters left behind

$digraph = array();
$inDegree = array();

foreach ($textAr as $line) {
    $line = preg_replace('#([0-9a-zA-Z_.\-]*[0-9a-zA-Z])[^0-9a-zA-Z]+([0-9a-zA-Z][0-9a-zA-Z_.\-]*)#', '$1 $2', $line);
    $arc = explode(" ", $line);
    echo '"'.$arc[0].'" -> "'.$arc[1].'";'."\n";
    if(!(array_key_exists($arc[0],$digraph))){
       $digraph[$arc[0]]=array();
    }
    $visited[$arc[0]]=0;
    $visited[$arc[1]]=0;
    if(!(array_key_exists($arc[1],$inDegree))){
       $inDegree[$arc[1]]=1;
    } else {
       $inDegree[$arc[1]]+=1;    
    }
    if(!(array_key_exists($arc[0],$inDegree))){
       $inDegree[$arc[0]]=0;
    }
    array_push($digraph[$arc[0]],$arc[1]);
}

$root = "r";
foreach($inDegree as $node => $degree){
    if($degree==0){
       $root=$node;
    }
}

function eNewick($source,$digraph,$inDegree,$internalNodesDisplay){
   global $reticulationVerticesFound,$visited,$reticulationVertices;
   $eNewick = "";

   // if $source is a reticulation vertex, compute its number
   if($inDegree[$source]>1){            
      if(!(array_key_exists($source,$reticulationVertices))){
         $reticulationVerticesFound += 1;
         $reticulationNumber = $reticulationVerticesFound;
         $reticulationVertices[$source] = $reticulationVerticesFound;
      } else {
         $reticulationNumber = $reticulationVertices[$source];
      }
   }
      
   if($visited[$source]==0){
      // if $source was not visited yet, recursively visit its children
      $visited[$source]=1;
      if(array_key_exists($source,$digraph)){
         $eNewick = "(";
         $i = 0;
         foreach($digraph[$source] as $index => $child){
            if($i > 0){
               $eNewick .= ",";
            }
            $eNewick .= eNewick($child,$digraph,$inDegree,$internalNodesDisplay);
            $i += 1;
         }         
         $eNewick .= ")";
      }
   }
   
   
   if(($internalNodesDisplay==1) or !(array_key_exists($source,$digraph))){
      $eNewick .= $source;
   }

   // if $source is a reticulation vertex, label it with its number
   if($inDegree[$source]>1){            
      $eNewick .= "#H".$reticulationNumber;
   }
   return $eNewick;
}

}

?>
}

</script>



<script src="js/viz.js"></script>


    <script>
      
      function inspect(s) {
        return "<pre>" + s.replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/\"/g, "&quot;") + "</pre>"
      }
      
      function src(id) {
        return document.getElementById(id).innerHTML;
      }
      
      function example(id, format, engine) {
        var result;
        try {
          result = Viz(src(id), format, engine);
          if (format === "svg")
            return result;
          else
            return inspect(result);
        } catch(e) {
          return inspect(e.toString());
        }
      }
      
      document.body.innerHTML += "<hr/><h1>Visualization of a rooted phylogenetic network</h1>";
      document.body.innerHTML += example("cluster", "svg");
      
    </script>
    

<hr/>

<?php 

$reticulationVerticesFound = 0;
$visited = array();
$reticulationVertices = array();
echo("<h2>eNewick string of the network</h2>

<p>You can use the codes below to visualize the phylogenetic network in <a href=\"http://dendroscope.org/\">Dendroscope</a> (menu <i>File</i>, <i>Enter Trees or Networks...</i>).</p>

<h3>With internal node names</h3><tt>".eNewick($root,$digraph,$inDegree,1).";</tt></p>");

$reticulationVerticesFound = 0;
$visited = array();
$reticulationVertices = array();
echo("<h3>Without internal node names</h3><tt>".eNewick($root,$digraph,$inDegree,0).";</tt></p>");

include("footer2.php");

?>


</body></html>