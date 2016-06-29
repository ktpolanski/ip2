### web service providing hypergeometric test for gene clusters against genes with motif ###

use Data::Dumper qw(Dumper);

# -------------------
# read in parameters
# -------------------

my $input_dir = $ARGV[0];
my $output_dir = $ARGV[1]; 

# -------------------
# Constants
# --------------------
#These should all be in $input_dir
my $cluster_file = 'Webapp_clusters.txt';
my $pvals_file = 'Webapp_p_values.txt';
my $which_genes_file = 'Webapp_which_genes.txt';
my $table_file = 'Webapp_table_gene_counts.txt';

my $logo_dir = $input_dir."/logos2";


my $test_for_large_overlaps = TRUE; #looking for over-representation

# --------------
# read in files
# --------------

#clusters, a hash of arrays. Cluster names are keys and arrays are list of gene ids in that cluster, use:
#@{$clusters{$all_cluster_IDs[$c]}}

# read in clusters into hash
my %clusters = ();
open(CLUSTERINPUTFILE, $input_dir.'/'.$cluster_file) or die "Cannot open '$cluster_file': $!";
while (<CLUSTERINPUTFILE>) {
 	chomp;
 	my ($cluster_id, $gene_id) = split("\t"); 
 	push @ { $clusters{$cluster_id} }, $gene_id;
}
close CLUSTERINPUTFILE;

#Read from p values file

#all_cluster_IDs, array of cluster names, use: $all_cluster_IDs[$c] 

#%num_genes_in_cluster, hash of values with keys being cluster names , use: $num_genes_in_cluster{$all_cluster_IDs[$c]}

#@p_values, motif * cluster array of p-values, used like $p_values[$motif][$cluster]

# pattern_model_IDs,  array of motif names, use: $pattern_model_IDs[$motif]


open(PVALSINPUTFILE, $input_dir.'/'.$pvals_file) or die "Cannot open '$pvals_file': $!";

#first line is cluster names
my $line = <PVALSINPUTFILE>;
$line =~ s/^\s+//;
my @all_cluster_IDs = split(/\s+/, $line);

#next line is sizes
my %num_genes_in_cluster;
$line = <PVALSINPUTFILE>;
$line =~ s/^\s+//;
my @sizes = split(/\s+/, $line);

for(my $c = 0; $c <= $#all_cluster_IDs; $c++){
	
	$num_genes_in_cluster{$all_cluster_IDs[$c]} = $sizes[$c];
	
}

#now motif names and p values
my @p_values = ();
my @pattern_model_IDs = ();

while ($line = <PVALSINPUTFILE>) {
 	chomp;
 	my @data = split(/\s+/, $line);
 	#first col is motif name
 	push @pattern_model_IDs, $data[0];
 	shift(@data);
 	#rest of line is pvals for that motif
 	push @p_values, \@data;
}
close PVALSINPUTFILE;

#read form which_genes file

my %gene_ids = ();

#gene_ids,  motif by cluster list of gene ids. Keys are motif, cluster. Values are array of gene ids.
#my %tmp_hash;
#@tmp_hash{@all_cluster_IDs} = undef
#@gene_ids{@pattern_model_IDs} = {%tmp_hash};
    
open(WHICHGENESFILE, $input_dir.'/'.$which_genes_file) or die "Cannot open '$which_genes_file': $!";
#ignore header and empty line

$line=<WHICHGENESFILE>;

while ($line = <WHICHGENESFILE>) {
 	chomp;
 	my @data = split(/\s+/, $line);
 	if ($#data>1){ #ignore empty line		
 		#first col is motif name, third is cluster id, 6th is comma seperated list of gene ids
 		#create array of gene ids
 		my $motif = $data[0];
 		my $cluster = $data[2];
 		my @gene_ids_for_motif_and_cluster = ();
 		if ($#data > 5){
 			@gene_ids_for_motif_and_cluster = split(/,/, $data[5]);	
 		}
 		$gene_ids{$motif}{$cluster} =  [ @gene_ids_for_motif_and_cluster]		
 		
 	}
}
 
close WHICHGENESFILE;

#read table counts file

#@universe_with_motif, list of scalars, key is motif name, vlaue number of each motif in universe,
#$size_of_universe, scalar, number of genes in universe

open(TABLEFILE, $input_dir.'/'.$table_file) or die ("Cannot open '$table_file': $!");
#ignore header
$line=<TABLEFILE>;
#next line starts with number of genes in universe
$line=<TABLEFILE>;
$line =~ s/^\s+//;
my @tmp = split(/\s+/, $line);
my $size_of_universe = $tmp[0];

my %universe_with_motif = ();

#following lines have motif name as col 1 and frequency in universe as col2
while ($line = <TABLEFILE>) {
 	chomp;
 	my @data = split(/\s+/, $line);
 	my $motif = $data[0];
 	my $num = $data[1];
	$universe_with_motif{$motif} = $num;	
}
 
close TABLEFILE;

##HTML output
if(! -d $output_dir){
	mkdir($output_dir, 0755)  or die("Error creating output directory");
}

mkdir($output_dir.'/which_genes', 0744) or die("Error creating output directory");

#make html files that pop up when clicking on p values. 
write_which_gene_ids_file($output_dir."/which_genes", \@p_values, \%gene_ids, \%universe_with_motif, $size_of_universe, \@all_cluster_IDs, \@pattern_model_IDs, \%num_genes_in_cluster, $test_for_large_overlaps);                   

my $htmlfile = $output_dir.'/hypergeometric_motif_test_result.html'; # full path to file

my $numclusters = @all_cluster_IDs;
my $numpssms = @pattern_model_IDs;

open (HTML, ">$htmlfile") or die "cannot write a results file at $htmlfile\n";

print HTML << "RESULTS" ;

<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">

<html>
<head>
<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.5/jquery.min.js"></script>
<script type="text/javascript">

var numClusters = $numclusters;
var numPSSMs = $numpssms;

<!-- code to show table with fixed column and row heading is modified from a jQuery script -->
<!-- downloaded from http://fixed-header-using-jquery.blogspot.com -->
<!-- This script adds the ability to add and remove rows and columns dynamically -->

function initialise(){
	
	fnAdjustTable();
	sizeTable();
	document.getElementById('tblContainer').style.visibility='visible';

}

var brow='mozilla';
	jQuery.each(jQuery.browser, function(i, val) {
		if(val==true){
			brow=i.toString();
		}
	});


fnAdjustTable=function(){
	

	var colCount=\$('#datarow0>td').length; //get total number of data columns
	var rowCount = \$('#firstcol tr').length;
	var firstVisibleRow=0;
	for(var n=0;n<rowCount;n++){
		var row = '#datarow'+n;
		if(\$(row).is(":visible")){
			firstVisibleRow	= n;
			break;
		}	
	}

	var visibleCol;
	for(m=0;m<colCount;m++){
		//sets width of col titles to match width of data cells
		var header = '#header'+m;
		if(\$(header).is(":visible")){
			//don't resize if col is hidden. This is important as when using Safari or Chrome,
			//reading the width of a td with display=none will result in
			//it being displayed again, just in the top row of the data table.
			
			//Size header cell using first visible data cell in that column
			var datacell = '#wm'+firstVisibleRow+'c'+m;
			\$(header).width(\$(datacell).width());			
			visibleCol = m; //note a visible col for later use
		}
	}
	
	
	for(var n=0;n<rowCount;n++){
		//sets height of row titles to match data cells
		var header = '#headerrow'+n;
		if(\$(header).is(":visible")){
			var datacell = '#wm'+n+'c'+visibleCol;
			var dataheight = \$(datacell).innerHeight();
			\$(header+' td').innerHeight(dataheight);
			\$(header+' img').attr('height', dataheight);
		}
	}
	
	
	
}


function sizeTable(){
	//resizes the outer table
	//initially called after fnAdjustTable
	
	var outerTable = document.getElementById("tblContainer");
	var tableLeft = outerTable.offsetParent.offsetLeft + outerTable.offsetLeft;
	var tableTop = outerTable.offsetParent.offsetTop + outerTable.offsetTop;
	
	//get window width, including any scroll bars as these will disappear after table resized
	if(typeof(window.innerWidth) == 'number'){
		//non-IE or IE9
		var pageWidth = window.innerWidth;
		var pageHeight = window.innerHeight;
	}else{
		//IE6-8
		var pageWidth = document.documentElement.clientWidth;
		var pageHeight = document.documentElement.clientHeight;	
	}

	var tableWidth = pageWidth-tableLeft-40;
	var tableHeight = pageHeight-tableTop-55;
	var firstColWidth = document.getElementById('firstTd').clientWidth;
	var firstRowHeight = document.getElementById('firstTd').clientHeight;

	
	var colHeaders = document.getElementById("divHeader");
	
	colHeaders.style.width = (tableWidth-firstColWidth)+"px";
	
	var rowHeaders = document.getElementById("firstcol");
	rowHeaders.style.height = (tableHeight-firstRowHeight)+"px";
	var dataTable = document.getElementById('table_div');
	
	//add extra size to allow for scrollbars
	dataTable.style.width = (16+tableWidth-firstColWidth)+"px";
	dataTable.style.height = (16+tableHeight-firstRowHeight)+"px";

	
	var leftcol = document.getElementById('leftcol');
	var height = (document.documentElement.clientHeight - leftcol.offsetTop-10) + "px";
	leftcol.style.height = height;

}


function syncScroll(){
	
	if(brow=='mozilla'){
		//Firefox won't reliably size the column headings and data table to the same widths
		//Data about 0.4px more per col
		document.getElementById("divHeader").scrollLeft = (document.getElementById("table_div").scrollLeft/document.getElementById("table_div").scrollWidth) * document.getElementById("divHeader").scrollWidth;
		document.getElementById("firstcol").scrollTop = (document.getElementById("table_div").scrollTop/document.getElementById("table_div").scrollHeight) * document.getElementById("firstcol").scrollHeight;
		
	}else{
		document.getElementById("divHeader").scrollLeft=document.getElementById("table_div").scrollLeft;
		document.getElementById("firstcol").scrollTop=document.getElementById("table_div").scrollTop;
	}
	
}


function rowHeaderScroll(){
	//called when row headings scrolled due to user using Find to search for a value that is found in the row headings
	//ensures that data table scrolled too.
	
	if(brow=='mozilla'){
		//Firefox won't reliably size the column headings and data table to the same widths
		//Data about 0.4px more per col
		document.getElementById("table_div").scrollTop = (document.getElementById("firstcol").scrollTop/document.getElementById("firstcol").scrollHeight) * document.getElementById("table_div").scrollHeight;
		
	}else{
		document.getElementById("table_div").scrollTop=document.getElementById("firstcol").scrollTop;
	}

	
}


function colHeaderScroll(){
	
	//as above but for when col headers scrolled
	
	if(brow=='mozilla'){
		//Firefox won't reliably size the column headings and data table to the same widths
		//Data about 0.4px more per col
		document.getElementById("table_div").scrollLeft = (document.getElementById("divHeader").scrollLeft/document.getElementById("divHeader").scrollWidth) * document.getElementById("table_div").scrollWidth;
		
	}else{
		document.getElementById("table_div").scrollLeft=document.getElementById("divHeader").scrollLeft;
	}

}


function filter()
{
	
		
	var pVal1 = document.getElementById('pVal1').value;
	var pVal2 = document.getElementById('pVal2').value;
	var pVal3 = document.getElementById('pVal3').value;

	
	//validate p values entered
	
	pVal1 = parseFloat(pVal1);
	if(isNaN(pVal1) || (pVal1 < 0 ) || (pVal1 > 1)){
		document.getElementById('pVal1').focus();
		alert('Please enter a valid probability value, 0 <= p <= 1');
		return;
	}
	pVal2 = parseFloat(pVal2);
	if(isNaN(pVal2) || (pVal2 < 0 ) || (pVal2 > 1)){
		document.getElementById('pVal2').focus();
		alert('Please enter a valid probability value, 0 <= p <= 1');
		return;
	}
	pVal3 = parseFloat(pVal3);
	if(isNaN(pVal3) || (pVal3 < 0 ) || (pVal3 > 1)){
		document.getElementById('pVal3').focus();
		alert('Please enter a valid probability value, 0 <= p <= 1');
		return;
	}
	
	if((pVal1 <= pVal2) || (pVal2 <= pVal3)){
		document.getElementById('pVal3').focus();
		alert('Probability thresholds should be increasingly less stringent (higher) in value ');
		return;
	}
	


	document.getElementById('progbox').style.display='block';
	document.body.style.cursor='wait';
	
	setTimeout("continueFilter("+pVal1+", "+pVal2+", "+pVal3+")", 1000);
	
	//delay allows display of progress box and cursor to be changed

	
}


function continueFilter(pVal1, pVal2, pVal3){
	
	var cellsInCol = new Array(numClusters);
	var cellsInRow = new Array(numPSSMs);
	var showRow;
	var numInCol;
	var visibleRows = 0;
	var visibleCols = 0;
	var numVisible = new Array(numClusters);

	
	var pColour1 = document.getElementById('pColour1').style.backgroundColor;
	var pColour2 = document.getElementById('pColour2').style.backgroundColor;
	var pColour3 = document.getElementById('pColour3').style.backgroundColor;
	
	var fontColours = new Array(3);
	fontColours[0] = (colorBrightness(pColour1) > 256)? "black":"white";
	fontColours[1] = (colorBrightness(pColour2) > 256)? "black":"white";
	fontColours[2] = (colorBrightness(pColour3) > 256)? "black":"white";
		
	for (var c = 0; c < numClusters; c++){
		cellsInCol[c] = 0;
	}
	
	
	//update matix
	for(var w = 0; w < numPSSMs; w++)
	{
		cellsInRow[w] = 0;
		for (var c = 0; c < numClusters; c++){
			var ctrl = document.getElementById('wm'+w+'c'+c);
			var p = parseFloat(ctrl.innerHTML);
			
			if(p > pVal1){
				ctrl.style.visibility = "hidden";
			}else{
				ctrl.style.visibility = "visible";
				cellsInCol[c]++;
				cellsInRow[w]++;
				if(p <= pVal3){
					ctrl.style.backgroundColor=pColour3;
					ctrl.style.color=fontColours[2];
				}else if (p <= pVal2){
					ctrl.style.backgroundColor=pColour2;
					ctrl.style.color=fontColours[1];
				}else{
					ctrl.style.backgroundColor=pColour1;
					ctrl.style.color=fontColours[0];
				}
			}	
		}
	}
	
	//hide rows with no values in
	for(var w = 0; w < numPSSMs; w++)
	{
		if(cellsInRow[w]){
			document.getElementById('headerrow' + w).style.display="";	//default is to show
			document.getElementById('datarow' + w).style.display="";	//default is to show
			visibleRows++;
		}else{
			document.getElementById('headerrow' + w).style.display="none";	//collapse space
			document.getElementById('datarow' + w).style.display="none";
		}	
	}
	
	//collapse cols with no entries. Only Firefox supports visible="collapse" for <col/>,
	//so have to use display on individual cells instead. Set to "none" will collapse them.
	//Cant do this on some cells in col, as some rows will be pushed right, other not,
	//destroying data integegrity. So just do this to cols were no cells needed.
	
	for (var c = 0; c < numClusters; c++){
		var ctrls = \$('td[name=col'+c+']');
		for (var cell = 0; cell < ctrls.length; cell++){
			if(cellsInCol[c]){
				ctrls[cell].style.display="";
			}else{
				ctrls[cell].style.display="none";
			}	
		}
		if(cellsInCol[c]){
			visibleCols++;
		}
	}
	
	//eliminating a row could shrink the width of a column if it contains the widest value in the column
	//Proabaly won't happen as all p values printed to same number of decimal places, but should be prepared
	//if it does. Conversely, putting a row back in could increase width.
	//Re-calculate column widths
	fnAdjustTable();
	sizeTable();
	
	//Same thing could happen when removing/adding a column if it contains the values in a row with thre most number
	//of lines, it will increase the height of the row. Shouldn't happen in this table as nothing 
	//should be split over more than one row, but including code to deal with this makes table more flexible
	
	document.getElementById('title').innerHTML = "Showing "+visibleRows+"/"+numPSSMs+" weight matrices and "+visibleCols+"/"+numClusters+" clusters";
	if(visibleRows == 0 || visibleCols == 0){
		document.getElementById('warning').style.display = 'block';
	}else{
		document.getElementById('warning').style.display = 'none';	
	}
	
	document.body.style.cursor='default';
	document.getElementById('progbox').style.display='none';
}

function colorBrightness(rgb){
	
	var digits = /rgb\\((\\d+),\\s*(\\d+),\\s*(\\d+)\\)/.exec(rgb);
	var r = parseInt(digits[1]);
	var g = parseInt(digits[2]);
	var b = parseInt(digits[3]);
	return (r+g+b);	
	
}


function clearfilter(){
	
	
	document.getElementById('pVal1').value = "";
	document.getElementById('pVal2').value = "";
	document.getElementById('pVal3').value = "";
	
	document.getElementById('progbox').style.display='block';
	document.body.style.cursor='wait';
	
	setTimeout("continueClearFilter()", 1000);
	
}

function continueClearFilter(){
	
	//show all rows and cols first
	for (var c = 0; c < numClusters; c++){
		var ctrls = \$('td[name=col'+c+']');
		for (var cell = 0; cell < ctrls.length; cell++){
			ctrls[cell].style.display="";
		}
	}

	for(var w = 0; w < numPSSMs; w++){
		document.getElementById('headerrow' + w).style.display="";	//default is to show
		document.getElementById('datarow' + w).style.display="";	
		for (var c = 0; c < numClusters; c++){
			var ctrl = document.getElementById('wm'+w+'c'+c);
			ctrl.style.visibility = "visible";
			ctrl.style.backgroundColor="white";
			ctrl.style.color="black";
		}
	}

	//re-adjust column widths
	fnAdjustTable();
	sizeTable();
	
	document.getElementById('warning').style.display = 'none';	
	document.getElementById('title').innerHTML = "Showing all "+ numPSSMs + " weight matrices and all " + numClusters + " clusters";
	document.body.style.cursor='default';
	document.getElementById('progbox').style.display='none';

}

</script>

<style type = "text/css">
			
	.headerTable {margin:0px;
					padding:0px;
					border-width:0px;
					border-style:none;}
	
	.headerTable td {padding:5px;
					margin:0px;
					font-size:14px;
					white-space:nowrap;
					overflow:hidden;
					text-align:center;
					background-color:powderblue;
					border:1px solid white;}
					
					
	.rowHeaderTable {margin:0px;
					padding:0px;
					border-width:0px;
					border-style:none;}

					
	.rowHeaderTable td {padding:5px;
					margin:0px;
					font-size:14px;
					white-space:nowrap;
					overflow:hidden;
					text-align:center;
					border:1px solid white;
					background-color:powderblue;}
					
			
	.dataTable {margin:0px;
					padding:0px;
					border-width:0px;
					border-style:none;}
	
	.dataTable td {padding:10px 5px 10px 5px;
					margin:0px;
					font-size:14px;
					white-space:nowrap;
					text-align:right;
					background-color:white;
					border:1px solid white;}
					
	.progbox {
			position:fixed;
			border:2px solid black;
			font-size:16px;
			font-weight:bold;
			width:240px;
			height:50px;
			top:50%;
			left:50%;
			margin-top:-25px;
			margin-left:-85px;
			background-color:white;
			cursor:wait;
			z-index:100}
					
	.colourpicker {
			position:fixed;
			border-width:1px;
			border-style:ridge;
			border-color:black;
			background-color:white;
			z-index:100;
			visibility:hidden;}
			
	.colourpicker table {
			border:0px none white;
			table-layout:fixed;
	}
	
				
	.colourpicker td {
			margin:0px;
			padding:0px;
			width:32px;
			height:20px;
			border-width:1px;
			border-style:solid;
			border-color:white;	
	}
	
	.colourpicker td:hover {
			border-color:#000000; 	
	}
	
	.detailsbox{
		background-color:white;
		border: 2px solid #287880;
		padding:0px;
		position:absolute;
		box-shadow:5px 5px 5px 2px #888888;		
	}
	
	.detailsbox table{
		border-width:0px;
		border-style:none;
		font-size:12px;
		margin:10px 0px 5px 10px;
	}
	
	.detailsbox td {
		padding:5px 0px 5px 5px;
		color:#287880;
		vertical-align:top;
	}
	
	.detailsbox img {
		border-width:0px;
		border-style:none;	
	}
				
				
</style>

</head>

<body style='font-family:sans-serif;background-color:white' onload='initialise()'>

<div id = 'progbox' class='progbox' style='display:block'>
	<div id='progtext' style='padding:15px'>Loading data ....</div>
</div>


<div id='div_titles'  style='margin-left:auto;margin-right:auto;padding-left:10px;padding-right:10px;font-size:16px;' align='center'><!-- showns input params -->
	<p align='center' style='font-family:sans-serif;font-size:24px;font-weight:bold'>Hypergeometric Motif Test results</p>

</div>


<div id = 'leftcol' style = 'position:fixed;left:0px;top:120px;width:250px;padding-left:10px;padding-right:10px;overflow-y:auto;font-size:14px;font-family:sans-serif;'><!-- left col -->

	<p><b>Filter values</b></p>
	<p><i>Filter and colour-code p-values and visualise the most important by setting the colours and thresholds below</i></p>
	<ul style='padding-left:15px'>
		<li>Set the first p-threshold to the most stringent (lowest) value. Any p-values &lt= to this will be shown in the selected colour</li>
		<li>Set the other thresholds to increasingly less stringent (higher) values</li>
		<li>The left column displays the logos repesented by the corresponding weight matrix, and the top row displays the cluster names plus the number of genes each contains</li>
		<li>Any p-values greater than the highest threshold will not be displayed</li>
		<li>Any weight matrix with no p-value displayed for any cluster will be hidden, eliminating that row</li>
		<li>Likewise, any cluster with no p-value displayed for any weight matrix will be hidden, eliminating that column</li>
		<li>Note that the p-values have not been corrected for multiple testing</li>
		<li>Hover the mouse pointer over a cell to display its contents below the table, or click on row or column headings for more information</li>
	</ul>
	
	<div style = 'position:absolute;top:560px;left:0px;width:60px;text-align:right;'>
		p &lt=
	</div>
	<div style = 'position:absolute;top:555px;left:65px;'>
		<input type='text' style='text-align:right;width:40px;' id='pVal3'/>
	</div>
	<div style="position:absolute;left:115px;cursor:pointer;top:555px;width:30px;height:20px;background-color:rgb(0,0,0);border:solid black 1px" id="pColour3" onclick="showColourPicker(this);" onmouseout="hideColourPicker(event)">
		
	</div>
	
	<div style = 'position:absolute;top:590px;left:0px;width:60px;text-align:right;'>
		p &lt=
	</div>
	<div style = 'position:absolute;top:585px;left:65px'>
		<input type='text' style='text-align:right;width:40px;' id='pVal2'/>
	</div>
	<div style="position:absolute;left:115px;cursor:pointer;top:585px;width:30px;height:20px;background-color:rgb(128,128,128);border:solid black 1px" id="pColour2" onclick="showColourPicker(this);" onmouseout="hideColourPicker(event)">
		
	</div>
	
	<div style = 'position:absolute;top:620px;left:0px;width:60px;text-align:right;'>
		p &lt=
	</div>
	<div style = 'position:absolute;top:615px;left:65px'>
		<input type='text' style='text-align:right;width:40px;' id='pVal1'/>
	</div>
	<div style="position:absolute;left:115px;cursor:pointer;top:615px;width:30px;height:20px;background-color:rgb(224,224,224);border:solid black 1px" id="pColour1" onclick="showColourPicker(this);" onmouseout="hideColourPicker(event)">
		
	</div>
	
	
	<div style = 'position:absolute;top:650px;left:0px;width:100px;'>
		<span style='position:absolute;right:5px'><button type='button' onclick='clearfilter()'>Clear</button></span>
	</div>
	<div style = 'position:absolute;top:650px;left:105px'>
		<button type='button' onclick='filter()'>Apply</button>
	</div>
	
	<!-- colour picker -->

	<div id = 'colourpicker' class='colourpicker' onmouseout='hideColourPicker(event)' onclick="selectColour(event, this)">
		<table>
			<tr>
				<td style='background-color:rgb(224,224,224)'> </td>
				<td style='background-color:rgb(255,128,128)'> </td>
				<td style='background-color:rgb(255,255,128)'> </td>
				<td style='background-color:rgb(128,255,128)'> </td>
				<td style='background-color:rgb(128,255,255)'> </td>
				<td style='background-color:rgb(128,128,255)'> </td>
				<td style='background-color:rgb(255,128,255)'> </td>
			</tr>
			<tr>
				<td style='background-color:rgb(128,128,128)'> </td>
				<td style='background-color:rgb(255,0,0)'> </td>
				<td style='background-color:rgb(255,255,0)'> </td>
				<td style='background-color:rgb(0,255,0)'> </td>
				<td style='background-color:rgb(0,255,255)'> </td>
				<td style='background-color:rgb(0,0,255)'> </td>
				<td style='background-color:rgb(255,0,255)'> </td>
			</tr>
			
			<tr>
				<td style='background-color:rgb(64,64,64)'> </td>
				<td style='background-color:rgb(128,0,0)'> </td>
				<td style='background-color:rgb(255,128,0)'> </td>
				<td style='background-color:rgb(0,128,0)'> </td>
				<td style='background-color:rgb(0,128,128)'> </td>
				<td style='background-color:rgb(0,0,128)'> </td>
				<td style='background-color:rgb(128,0,128)'> </td>
			</tr>
			<tr>
				<td style='background-color:rgb(0,0,0)'> </td>
				<td style='background-color:rgb(64,0,0)'> </td>
				<td style='background-color:rgb(128,64,0)'> </td>
				<td style='background-color:rgb(0,64,0)'> </td>
				<td style='background-color:rgb(0,64,64)'> </td>
				<td style='background-color:rgb(0,0,64)'> </td>
				<td style='background-color:rgb(64,0,64)'> </td>
			</tr>
			
		</table>
	</div>

</div>

<div id = 'div_data' style = 'position:absolute;left:270px;top:70px;padding:10px;font-family:sans-serif;'><!-- p=vals -->

	<div style="float:right;padding:5px;margin-right:25px">
		<button id="zoomout" type="button" onclick="zoom(-3)">Zoom out</button><button id="zoomin" type="button" onclick="zoom(3)">Zoom in</button>
	</div>
	<div id = 'div_display_info' style = 'float:left;padding-left:10px;padding-right:10px;font-size:14px;font-family:sans-serif;'>
		<p>
			<span id = 'title' style='font-style:italic'>Showing all $numpssms weight matrices and all $numclusters clusters</span>
		</p>
	</div><br/>

	<span  id = 'warning' style='text-align:center;font-size:20px;color:red;display:none'>No p-values exceed the probability threshold</span>

	<table id='tblContainer' style='visibility:hidden' cellspacing='0' cellpadding='0' border='0'>
		<tr>
			<td id='firstTd' valign='top'><button type='button' onclick='fullScreen()'><div id='fullsizearrow' style='font-size:18px;font-weight:900;transform:rotate(45deg);-ms-transform:rotate(45deg);-moz-transform:rotate(45deg);-webkit-transform:rotate(45deg);-o-transform:rotate(45deg);'>&#8592</div></button></td>
			<td valign='bottom'>
				<div id='divHeader' style='overflow:hidden;'  onscroll='colHeaderScroll()'>
					<table class = 'headerTable' cellspacing='0' cellpadding='0'>
						<tr>
						
RESULTS

#print col headings (= cluster names)
for (my $c=0;$c<=$#all_cluster_IDs;$c++) {
	print HTML "\t\t\t\t\t\t\t<td onclick = 'showClusterDetails(this)' name = 'col".$c."' id = 'col".$c."' onmouseover='highlightTitles(-1,".$c.");' onmouseout='unHighlightTitles(event, -1,".$c.");'  style='cursor:pointer'><div id='header".$c."'>".$all_cluster_IDs[$c]."<br/>".$num_genes_in_cluster{$all_cluster_IDs[$c]}." genes</div></td>\n";

	print HTML "\t\t\t\t\t\t<script type = 'text/javascript'>\n";
	print HTML "\t\t\t\t\t\t\tvar geneList = new Array();\n";
	print HTML "\t\t\t\t\t\t\tgeneList.push('".$all_cluster_IDs[$c]."');\n";
	
	my @gene_list = @{$clusters{$all_cluster_IDs[$c]}};
	for(my $g=0; $g<=$#gene_list;$g++){
		print HTML "\t\t\t\t\t\t\tgeneList.push('".$gene_list[$g]."');\n";
	}
	
	print HTML "\t\t\t\t\t\t\t\$('#col".$c."').data('genes', geneList);\n";
	print HTML "\t\t\t\t\t\t</script>\n";

}
																								
print HTML << "RESULTS" ;

						</tr>
					</table>
				</div>
			</td>
		</tr>
		<tr>
			<td valign='top'>
				<div id='firstcol' style='overflow:hidden;' onscroll='rowHeaderScroll()'>
					<table cellspacing='0' cellpadding='0' class = 'rowHeaderTable'>
					
RESULTS

#print row headings, = pssm names
system("mv",  $logo_dir, $output_dir."/img");

for (my $r=0;$r<=$#pattern_model_IDs;$r++) {
	
	my $motif_name = $pattern_model_IDs[$r];
	my $logoname = 'img/'.$motif_name.".png";
			
	print HTML "\t\t\t\t\t\t<tr id = 'headerrow".$r."' onclick='showPSSMDetails(this)' onmouseover='highlightTitles(".$r.",-1);' onmouseout='unHighlightTitles(event, ".$r.",-1);' style='cursor:pointer'><td id = 'firstcol".$r."'>".$motif_name."</td><td style='padding:0px'><img alt = '".$motif_name."' border='0'  name = '".$motif_name."' src = '".$logoname."' /></td></tr>\n";

}
					
print HTML << "RESULTS" ;
						
					</table>
				</div>
			</td>
			<td valign='top'>
				<div id='table_div' style='overflow:scroll;position:relative' onscroll='syncScroll()'>
					<table class = 'dataTable' cellspacing='0' cellpadding='0'>
					
RESULTS

#print p values
my $w = 0;
foreach my $array_ref (@p_values) {	#for each wm

    my @one_p_value_array = @{$array_ref};	#a row of p values
    print HTML "\t\t\t\t\t\t<tr id = 'datarow".$w."'>";
    for (my $c=0;$c<=$#one_p_value_array;$c++) {
		print HTML "<td id='wm".$w."c".$c."' name='col".$c."' onmouseover='highlightTitles(".$w.",".$c.");' onmouseout='unHighlightTitles(event, ".$w.",".$c.");' onclick = 'showWhichGeneIDs(event, ".$w.", ".$c.");' style='cursor:pointer'>".$one_p_value_array[$c]."</td>";
    
    
    }
    print HTML "</tr>\n";
    $w++;
}
#
#system("cp", "-R" , $mydir.'/img', $output_dir);#this is just arrow images for main webpage

						
print HTML << "RESULTS" ;						
						
					</table>
				</div>
			</td>
		</tr>
	</table>		
	<!--/div-->
	<div id='status_bar' style='font-size:14px;font-family:sans-serif;font-style:italic;float:right;margin:5px 25px 5px 0px'>
		
	</div>
</div>


<div id='pssm_details' class='detailsbox' style='display:none' onmouseout='hidePSSMDetails(event);'>
	<div style='width:100%;background-color:#287880;color:white;overflow:auto;'>
		<span id='technical_name' style='padding-left:2px'></span>
	</div>
	<div style='float:left'>
		<img id='logo_image' />
	</div>
</div>

<div id='cluster_details' class='detailsbox' style='display:none;overflow:hidden;' onmouseout='hideClusterDetails(event);' >
	<div style='width:100%;background-color:#287880;color:white;overflow:auto;'>
		<span id='cluster_name' style='padding-left:2px'></span>
	</div>
	<div  id = 'cluster_title' style='font-weight:bold;font-size:12px;margin:10px 0px 5px 10px;padding:5px 0px 5px 5px;color:#287880;'></div>
	<div id = 'cluster_genes' style='overflow-y:scroll;padding-bottom:5px'>
		<table border='0' id = 'cluster_info_table' style='padding-right:10px;'>
		
	
		</table>
	</div>

</div>


<script type='text/javascript'>

	function zoom(direction){
		
		var size = \$("#wm0c0" ).css("font-size");
		size = parseInt(size);
		
		var scrollPos = document.getElementById("table_div").scrollLeft/document.getElementById("table_div").scrollWidth;
		size = (size+direction);
		if ((size >= 2) && (size <= 20)){
			\$("#tblContainer td" ).css("font-size", size+"px");	
		}
		if(size > 2){
			document.getElementById('zoomout').disabled = false;
		}else{
			document.getElementById('zoomout').disabled = true;
		}
		if(size  < 20){
			document.getElementById('zoomin').disabled = false;
		}else{
			document.getElementById('zoomin').disabled = true;
		}
		
		fnAdjustTable();
		sizeTable();
		
		//adjust scroll, trying to key pa values in roughly the same place
		document.getElementById("table_div").scrollLeft = scrollPos*document.getElementById("table_div").scrollWidth;
		
		syncScroll();

	}

	document.getElementById('progbox').style.display='none';
	document.getElementById('progtext').innerHTML = "Updating ...";
	
	window.onresize = function() {sizeTable();}
	
	var colourToSet = null;
		
	function highlightTitles(w, c){
		
		var clusterName = "";
		var wmName = "";
		var pVal = "";
		
		if(c>=0){
			document.getElementById('col'+c).style.backgroundColor='#287880';
			document.getElementById('col'+c).style.color='white';
			
			var name = document.getElementById('header'+c).innerHTML;
			name = name.substr(0, name.indexOf("<br"));
			clusterName = "Cluster " + name;

			if(brow == 'mozilla'){
				var correction = document.getElementById('wm0c'+c).offsetLeft-document.getElementById('col'+c).offsetLeft;
				document.getElementById("divHeader").scrollLeft = document.getElementById("table_div").scrollLeft-correction;
			}
			
		}
		if(w>=0){
			\$("#headerrow"+w+" td:eq(0)").css('background-color', '#287880');
			\$("#headerrow"+w+" td").css('color', 'white');
			
			wmName =  "Weight Matrix " + \$('#headerrow'+w+' img').first().attr('name');
			if(brow == 'mozilla'){
				var correction  = document.getElementById('wm'+w+'c0').offsetTop-document.getElementById('firstcol'+w).offsetTop;
				document.getElementById("firstcol").scrollTop = document.getElementById("table_div").scrollTop-correction;
			}

		}
		if(c>=0 && w>=0){
			document.getElementById('wm'+w+'c'+c).style.borderColor='black';
			pVal = " : "+document.getElementById('wm'+w+'c'+c).innerHTML;
			clusterName +=", ";
		}
				
		document.getElementById('status_bar').innerHTML = (clusterName+wmName+pVal);

	}
	
	function unHighlightTitles(e, w, c){
		
		if(c>=0){
			
			var clearColTitle = true;
			
			if(\$('#cluster_details').is(":visible")){
			
				if (!e) {
					var e = window.event;
				}
			
				var toElement = (e.relatedTarget)? e.relatedTarget : e.toElement;
				while (toElement.nodeName != "BODY"){
					if (toElement == document.getElementById('cluster_details')){
						clearColTitle = false;
						break;
					}	
					toElement = toElement.parentNode;
				}
				if(clearColTitle == true){
					hideClusterDetails(e);
				}
			}
			if(clearColTitle == true){
				document.getElementById('col'+c).style.color='black';
				document.getElementById('col'+c).style.backgroundColor='powderblue';
			}

		}
		if(w>=0){
			
			var clearRowTitle = true;
			
			if(\$('#pssm_details').is(":visible")){
			
				if (!e) {
					var e = window.event;
				}
			
				var toElement = (e.relatedTarget)? e.relatedTarget : e.toElement;
				while (toElement.nodeName != "BODY"){
					if (toElement == document.getElementById('pssm_details')){
						clearRowTitle = false;
						break;
					}	
					toElement = toElement.parentNode;
				}
				if(clearRowTitle == true){
					hidePSSMDetails(e);
				}
			}
			if(clearRowTitle == true){
				\$("#headerrow"+w+" td:eq(0)").css('background-color', 'powderblue');
				\$("#headerrow"+w+" td:eq(1)").css('background-color', 'powderblue');
				\$("#headerrow"+w+" td").css('color', 'black');
			}
		}
		if(c>=0 && w>=0){
			document.getElementById('wm'+w+'c'+c).style.borderColor='white';
		}
		
		document.getElementById('status_bar').innerHTML = "";
		
	}


	
	function showColourPicker(col){
			
		var leftpos = col.offsetLeft+col.offsetParent.offsetLeft-2;

		var toppos = col.offsetTop+col.offsetParent.offsetTop-col.offsetParent.scrollTop - document.getElementById('colourpicker').offsetHeight;
		
		document.getElementById('colourpicker').style.left = leftpos+"px";
		document.getElementById('colourpicker').style.top = toppos+"px";
		document.getElementById('colourpicker').style.visibility="visible";
		
		colourToSet = col;
		
	}
	
	function hideColourPicker(e){
		
		
		if (!e) {
			var e = window.event;
		}
		var sourceElement = (e.srcElement)? e.srcElement:e.target;
		
		//first case, mouse leaves picker
		while(sourceElement.nodeName != "BODY"){
			if(sourceElement == document.getElementById('colourpicker')){
				//mouse moved from picker, or one of its children
				var toElement = (e.relatedTarget)? e.relatedTarget : e.toElement;
				while (toElement.nodeName != "BODY"){
					if (toElement == sourceElement){
						//mouse moved within picker, ignore
						return;	
					}	
					toElement = toElement.parentNode;
				}
				//if we get here, mouse moved from within picker to outside pickder
				colourToSet = null;
				document.getElementById('colourpicker').style.visibility="hidden";
				return;	
			}
			sourceElement = sourceElement.parentNode;
		}		
		//second case, mouse leaves colour box that activated picker
		//hide picker unless mouse moved to picker
		var toElement = (e.relatedTarget)? e.relatedTarget : e.toElement;
		while (toElement.nodeName != "BODY"){
			if (toElement == document.getElementById('colourpicker')){
				//mouse moved to picker, ignore
				return;	
			}	
			toElement = toElement.parentNode;
		}
		//if we get here, mouse moved from colour box, but not to picker
		colourToSet = null;
		document.getElementById('colourpicker').style.visibility="hidden";
		return;	
		
	}
	
	function selectColour(e,col){
			
		if (colourToSet){
			
			if (!e) {
				var e = window.event;
			}
			var sourceElement = (e.srcElement)? e.srcElement:e.target;
					
			if(sourceElement.nodeName == "TD"){
				colourToSet.style.backgroundColor = sourceElement.style.backgroundColor;	
				colourToSet = null;
			}
		}
		
		col.style.visibility="hidden";	
		
	}
	
	
	function fullScreen(){

		if(\$('#div_titles').is(":visible")){
			document.getElementById('div_titles').style.display='none';
			document.getElementById('leftcol').style.display='none';
			document.getElementById('div_data').style.left="10px";
			document.getElementById('div_data').style.top="10px";
			\$("#fullsizearrow").css('transform', 'rotate(225deg)');
			\$("#fullsizearrow").css('-ms-transform', 'rotate(225deg)');
			\$("#fullsizearrow").css('-moz-transform', 'rotate(225deg)');
			\$("#fullsizearrow").css('-webkit-transform', 'rotate(225deg)');
			\$("#fullsizearrow").css('-o-transform', 'rotate(225deg)');
			
		}else{
			document.getElementById('div_titles').style.display='block';
			document.getElementById('leftcol').style.display='block';
			document.getElementById('div_data').style.left="270px";
			document.getElementById('div_data').style.top="70px";
			\$("#fullsizearrow").css('transform', 'rotate(45deg)');
			\$("#fullsizearrow").css('-ms-transform', 'rotate(45deg)');
			\$("#fullsizearrow").css('-moz-transform', 'rotate(45deg)');
			\$("#fullsizearrow").css('-webkit-transform', 'rotate(45deg)');
			\$("#fullsizearrow").css('-o-transform', 'rotate(45deg)');

		}
		fnAdjustTable();
		sizeTable()	
		
	}
	
	function showPSSMDetails(currentrow){
		
		//find row 
		var row = \$(currentrow);
		
		var logoname = row.find("img").first().attr('src');	

		\$('#logo_image').attr('src', logoname);	
	
		\$('#pssm_details').data('row', row.attr('id'));
		
		var height = \$(row).outerHeight();
		var pos = \$(row).children("td").eq(1).offset();
		var left = pos.left;
		var top = pos.top;
		\$('#pssm_details').css('left', left+'px');
		
		pos = \$('#tblContainer').offset();
		var rightEdge = pos.left + \$('#tblContainer').outerWidth();
		var maxWidth = rightEdge-left;
		\$('#pssm_details').css('max-width', maxWidth+'px');
		
		pos =  \$('#firstcol').offset();
		var bottom = pos.top + \$('#firstcol').outerHeight();
		
		\$('#pssm_details').css('visibility', 'hidden');
		\$('#pssm_details').css('display', 'block');
		var imgW = \$('#logo_image').outerWidth();
		
		var newImgW = Math.min(maxWidth-10, 400);
		\$('#logo_image').css('width', newImgW +'px');	
		\$('#logo_image').css('height', (newImgW/2) +'px');
		var boxHeight = \$('#pssm_details').outerHeight();
		\$('#pssm_details').css('display', 'none');
		\$('#pssm_details').css('visibility', 'visible');
		
		top = Math.min(top, (bottom-boxHeight));
		\$('#pssm_details').css('top', top+'px');
		
		\$('#pssm_details').animate({width: 'show'}, 500);
		
	}
	
	function hidePSSMDetails(e){
	
		//mouse leaves row header or details box
		
		if(\$('#pssm_details').is(":not(:visible)")){
			return;
		}
		var row = \$('#pssm_details').data('row');
		
		if (!e) {
			var e = window.event;
		}
		var toElement = (e.relatedTarget)? e.relatedTarget : e.toElement;
		while (toElement.nodeName != "BODY"){
			if (toElement == document.getElementById('pssm_details') || toElement == document.getElementById(row)){
				//mouse moved within box, or back to row that activated it, so ignore
				return;	
			}	
			toElement = toElement.parentNode;
		}
		
		
		//if we get here, mouse moved from within box to outside box or row
		\$("#"+row+" td:eq(0)").css('background-color', 'powderblue');
		\$("#"+row+" td").css('color', 'black');
		\$('#pssm_details').data('row', null);
		\$("#pssm_details").animate({width: "hide"}, 500);
		
	}
	
	function showClusterDetails(currentCol){
		
		var genes = \$(currentCol).data('genes');
		var name = \$(currentCol).html();

		var clusterName = genes[0];
		var n = genes.length-1;
		
		\$("#cluster_name").html(clusterName);
		\$("#cluster_title").html("The cluster "+clusterName+" contains the following "+n+" genes:");
		
		while(document.getElementById('cluster_info_table').rows.length>0){
			document.getElementById('cluster_info_table').deleteRow(-1);
		}

		for(var g=1; g<genes.length; g++){
			if(!((g-1)%4)){
				var newrow = document.getElementById('cluster_info_table').insertRow(-1);
				newrow.insertCell(0);
				newrow.insertCell(1);
				newrow.insertCell(2);
				newrow.insertCell(3);
			}
			var row = document.getElementById('cluster_info_table').rows[document.getElementById('cluster_info_table').rows.length-1];
			row.cells[(g-1)%4].innerHTML = genes[g];
		}

		var width = Math.min(50, \$(currentCol).outerWidth());
		var pos = \$(currentCol).offset();
		var left = pos.left + (\$(currentCol).outerWidth()-width)/2;
		var top = pos.top + \$(currentCol).outerHeight();
		\$('#cluster_details').css('top', top+'px');
		
		pos = \$('#tblContainer').offset();
		var rightEdge = pos.left + \$('#tblContainer').outerWidth();
		left = Math.min(left, (rightEdge-\$('#cluster_details').outerWidth()));
				
		\$('#cluster_details').css('left', left+'px');
		
		var bottomEdge = pos.top + \$('#tblContainer').outerHeight();
		var maxHeight = bottomEdge-top;
		\$('#cluster_details').css('max-height', maxHeight+'px');
		\$('#cluster_genes').css('max-height', (maxHeight-60)+'px');

		\$('#cluster_details').data('col', currentCol.id);
		\$('#cluster_details').animate({height: 'show'}, 500);
		
	}
	
	function hideClusterDetails(e){
	
		//mouse leaves row header or details box
		
		if(\$('#cluster_details').is(":not(:visible)")){
			return;
		}
		var col = \$('#cluster_details').data('col');
		
		if (!e) {
			var e = window.event;
		}
		var toElement = (e.relatedTarget)? e.relatedTarget : e.toElement;
		while (toElement.nodeName != "BODY"){
			if (toElement == document.getElementById('cluster_details') || toElement == document.getElementById(col)){
				//mouse moved within box, or back to row that activated it, so ignore
				return;	
			}	
			toElement = toElement.parentNode;
		}
		
		
		//if we get here, mouse moved from within box to outside box or row
		\$("#"+col).css('background-color', 'powderblue');
		\$("#"+col).css('color', 'black');
		\$('#cluster_details').data('col', null);
		\$("#cluster_details").animate({height: "hide"}, 500);
		
	}
	
	function showWhichGeneIDs(e, motif, cluster){

    	if (!e) {
        	var e = window.event;
        }

        var top = e.screenY;
        var left = e.screenX-600;

        try{
             window.open("which_genes/" + motif + "_" + cluster + ".html", "", "top="+top+",left="+left+",menubar=no,status=no,toolbar=no,width=600,height=400,scrollbars=yes");
        }catch(err){
             window.open("which_genes/" + motif + "_" + cluster + ".html", "", "menubar=no,status=no,toolbar=no,width=600,height=400,scrollbars=yes");
    	}
    }


</script>
	
</body>

</html>
	
RESULTS

close(HTML);


sub write_which_gene_ids_file {

   #creates an html file that lists all the gene ids associated with a specific motif
   #One file created for each p-value, that this motif * cluster combination
   
   #write_which_gene_ids_file($output_dir."/which_genes", \@p_values, \%gene_ids, \%universe_with_motif, $size_of_universe, \@all_cluster_IDs, \@pattern_model_IDs, \%num_genes_in_cluster, $test_for_large_overlaps);                   

    my $file_path = $_[0];          	#directory to save files in
    my @p_values = @{$_[1]};        	#motif * cluster array of p-values, , ordered as in pvals input file
    my %gene_ids = %{$_[2]};        	#motif by cluster hash of gene id arrays, $gene_ids{$motif}{$cluster}
    my %universe_with_motif = %{$_[3]};	#hash of scalars, one for each motif
    my $size_of_universe = $_[4];       #number of genes in universe
    my @cluster_names = @{$_[5]};		#Names of the clusters, ordered as in pvals input file
    my @motif_names = @{$_[6]};			#Names of the motifs, ordered as in pvals input file
    my %size_of_cluster = %{$_[7]};		#Hash of cluster sizes. Names are keys
    my $test_for_over_representation = $_[8];
    
    #Files must be named m_c.html, where m is a zero-based row index of p_values representing motif, and
    #c is a zero based column index of p_values representing cluster. Motif and cluster name arrays are read from
    # p vals file so list these names in the corresponding order
        
     
    if (substr($file_path, -1, 1) ne "/"){
    	$file_path = $file_path."/";
    }
    
    for(my $motif = 0; $motif <= $#motif_names; $motif++){ #for each table row (motif)
    
    		
    		for(my $cluster = 0; $cluster <= $#cluster_names; $cluster++){ #for each cluster (column)
    		
    			my $m_name = $motif_names[$motif];
    			my $c_name = $cluster_names[$cluster];
    		
    			my $title = "Motif ".$m_name." Cluster ".$c_name;
				my @gene_ids_for_motif_and_cluster  = @{$gene_ids{$m_name}{$c_name}}; 
				my $num_genes_in_cluster_with_motif = $#gene_ids_for_motif_and_cluster+1;
    		
				my $percent_universe = sprintf("%.2g", $universe_with_motif{$m_name}/$size_of_universe);
				my $percent_cluster = sprintf("%.2g", $num_genes_in_cluster_with_motif/$size_of_cluster{$c_name});
			
    			my $htmlfile = $file_path.$motif."_".$cluster.".html"; #file for motif/cluster combination   

    		
    			open (FH, ">$htmlfile") or die "cannot write a results file at $htmlfile\n";
    			print FH << "RESULTS";
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">

<html>
<head>
        <title>$title</title>

        <script type="text/javascript">
        

                window.onresize = function(){
                        var height = (document.documentElement)? document.documentElement.clientHeight : document.clientHeight;
                        document.getElementById('scrollBox').style.height = (height - document.getElementById('rightCol').offsetTop-50)+"px";
                        var width =  (document.documentElement)? document.documentElement.clientWidth : document.clientwidth;
                        var tablewidth = (width - document.getElementById('rightCol').offsetLeft-20);
                        tablewidth = Math.min(tablewidth, 280);
                        document.getElementById('scrollBox').style.width = tablewidth + "px";
                }
                
                
        </script>

        <style type = "text/css">

                body    {
                                font-family:sans-serif;
                                background-color:white;
                                padding:10px;
                        }

                table.tbl   {
                                border-width:1px;
                                border-style:solid;
                                font-size:14px;
                                table-layout:fixed;
                                white-space:nowrap;
                        }

                .tbl th      {
                                padding:5px;
                                vertical-align:top;
                                font-weight:bold;
                                border-width:1px;
                                border-style:solid;
                }

                td      {
                                padding:5px;
                              	border-style:none;
                }
                
                .tbl td {
                	padding: 5px;
                	border-width:1px;
                	border-style:solid;
                	text-align:right;	
                }

                .leftCol {

                        left:10px;
                        top:100px;
                        position:fixed;
                        padding-left:10px;
                        width:320px;
                        height:205px;
 
                }

                .rightCol {
                        top:100px;
                        width:240px;
                        position:fixed;
                        left:345px;
                        height:205px;
   
                }

                .scrollBox {
                        overflow-y:auto;
                        overflow-x:hidden;
                        height:200px;
                }


        </style>

</head>

<body>
        <div style='font-size:18px;background-color:#CCD8DD;padding:10px;'>
            $title
        </div>
        <div style='font-size:16px;font-style:italic;margin-top:10px;'>
            The distribution of this motif is as follows
        </div>
        
        <div id='leftCol' class='leftCol'>
        	<table class = "tbl">
        		<tr>
        			<th></th><th>Universe</th><th>Cluster $cluster_names[$cluster]</th>
        		</tr>
        		<tr>
        			<th align="left">Size of cluster</th><td>$size_of_universe</td><td>$size_of_cluster{$c_name}</td>
        		</tr>
        		<tr>
        			<th align="left">Frequency of motif</th><td>$universe_with_motif{$m_name}</td><td>$num_genes_in_cluster_with_motif</td>
        		</tr>
        		<tr>
        			<th align="left">Proportion hits</th><td>$percent_universe</td><td>$percent_cluster</td>
        		</tr>
        	</table>
				
        	<div style='font-size:14px;margin-top:10px;'>
RESULTS
        	
        	if($test_for_over_representation){
        		print FH "\t\t\t\t<i>p</i>-value for over-representation in this cluster = ", sprintf("%.2g", $p_values[$motif][$cluster]), "\n";
        	}else{
        		print FH "\t\t\t\t<i>p</i>-value for under-representation in this cluster = ", sprintf("%.2g", $p_values[$motif][$cluster]), "\n";
        	}
        	
  			print FH << "RESULTS" ;

			</div>
		</div>
  		<div id='rightCol' class='rightCol'>
  			<div style='font-size:14px;margin-top:10px;margin-bottom:10px;font-weight:bold'>
				Cluster gene IDs with this motif
        	</div> 
        	<div class='scrollBox' id='scrollBox'>	
				<table style="border-style:none;font-size:14px">
        		
RESULTS
			for(my $gene_id = 0; $gene_id <= $#gene_ids_for_motif_and_cluster; $gene_id++){
				
				print FH "\t\t\t\t\t<tr><td>", $gene_ids_for_motif_and_cluster[$gene_id], "</td></tr>\n";
				
			}


        	print FH << "RESULTS" ;
        		
        		</table>
        	</div>
       </div>	
</body>
</html>

RESULTS


        	close(FH);
		} #end of cluster
    } #end of motif

} #end of function

