<?php

#$_GET['id']='76224_Vv';
# get info (annotation, v3 name etc.)
function annotation() {
  $handle = fopen('fantom-anno.csv', 'r');
  while(($line = fgetcsv($handle, 0, ',')) != FALSE) {
    $anno[$line[0]] = join('|', array_slice($line, 1));
  }
  echo "annotation loaded\n";
  fclose($handle);
  return $anno;
}

#
  array_shift($argv);
  $j = 0;
  foreach ($argv as $file) {
    if(file_exists($file)) {
      $i = 0; 
      if($j == 0) $anno = annotation();
      $content = file($file);
      foreach($content as $line) {
        if($i == 0) {
          $a = explode(' ', chop($line));
          $tcode = substr($a[3], 0, 7);
          list($iso) = explode('|', $anno[$tcode]);
          list($pcim_id) = explode('_', $file); 
          $handle = fopen($filename = $pcim_id."_Hs_$iso.csv", 'w');
          echo $filename ."\n";
          fwrite($handle, join(',', $a).',"'.$iso.'"'.PHP_EOL);
        }
        else if($i == 1) {
          fwrite($handle, 'rank,ID,Fabs,Frel,association_with_transcript,entrezgene_id,hgnc_id,uniprot_id,gene_name,description,type'.PHP_EOL);
        
        } else {
          $a = explode(',', $line);
          $a[1][0] = 'T';
          $a[3] = number_format($a[3], 4);
          array_pop($a); # remove ??
          fwrite($handle, join(',', $a) . ',');
          fputcsv($handle, explode('|', $anno[$a[1]]), ',', '"');
        }
        $i++;
      }
      fclose($handle);
    } else {
      echo "$file: not found\n";
    }
    $j++;
  }
?>
