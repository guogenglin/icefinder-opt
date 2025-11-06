/* Copyright 2012 Kai Blin. Licensed under the Apache License v2.0, see LICENSE file */

var svgene = {
    version: "0.1.6",
    label_height: 14,
    extra_label_width: 70,
    unique_id: 0
};

svgene.geneArrowPoints = function (orf, height, offset, border, scale) {
  var top_ = offset + svgene.label_height + border;
  var bottom = offset + svgene.label_height + height - border;
  var middle = offset + svgene.label_height + (height/2);
  if (orf.strand == 1) {
      var start = scale(orf.start);
      var box_end = Math.max(scale(orf.end) - (2*border), start);
      var point_end = scale(orf.end);
      points  = "" + start + "," + top_;
      points += " " + box_end + "," + top_;
      points += " " + point_end + "," + middle;
      points += " " + box_end + "," + bottom;
      points += " " + start + "," + bottom;
      return points;
  }
  if (orf.strand == -1) {
      var point_start = scale(orf.start);
      var end = scale(orf.end);
      var box_start = Math.min(scale(orf.start) + (2*border), end);
      points = "" + point_start + "," + middle;
      points += " " + box_start + "," + top_;
      points += " " + end + "," + top_;
      points += " " + end + "," + bottom;
      points += " " + box_start + "," + bottom;
      return points;
  }
};

svgene.ttaCodonPoints = function (codon, height, offset, border, scale) {
    var top_ = offset + svgene.label_height + height;;
    var bottom = offset + (2 * svgene.label_height) + height - border;
    var tip = Math.floor(scale(codon.start), scale(codon.end));
    var points = "" + tip + "," + top_;
    points += " " + (tip - 5) + "," + bottom;
    points += " " + (tip + 5) + "," + bottom;
    return points;
}

svgene.drawOrderedClusterOrfs = function(cluster, chart, all_orfs, borders, tta_codons,
                                         scale, i, idx, height, width,xAxis,
                                         single_cluster_height, offset) {
  chart.append("line")
    .attr("x1", 0)
    .attr("y1", (single_cluster_height * i) + svgene.label_height + (height/2))
    .attr("x2", width)
    .attr("y2", (single_cluster_height * i) + svgene.label_height + (height/2))
    .attr("class", "svgene-line");
	
	  chart.append('g')
    .attr('transform', 'translate('+0+',' + 100 + ')')
    .attr('class', 'axis')
    .call(xAxis);
	
	
  chart.selectAll("rect")
    .data(borders)
  .enter().append("rect")
    .attr("x", function(d){ return scale(d.start)})
    .attr("y", function(d){ var y = single_cluster_height * i  + offset; if (d.tool == "clusterfinder") y -= 3; return y})
    .attr("height", 5)
    .attr("width", function(d){ var width = scale(d.end) - scale(d.start); return width})
    .attr("class", function(d){ return "svgene-border-" + d.tool});
  chart.selectAll("polygon")
    .data(all_orfs)
  .enter().append("polygon")
    .attr("points", function(d) { return svgene.geneArrowPoints(d, height, (single_cluster_height * i), offset, scale); })
    .attr("class", function(d) { return "svgene-type-" + d.type + " svgene-orf"; })
    .attr("id", function(d) { return idx + "-cluster" + cluster.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-orf"; })
    .attr("style", function(d) { if (d.color !== undefined) { return "fill:" + d.color; } });
  chart.selectAll("polyline.svgene-tta-codon")
    .data(tta_codons)
  .enter().append("polyline")
    .attr("points", function(d) { return svgene.ttaCodonPoints(d, height, (single_cluster_height * i), offset, scale) })
    .attr("class", "svgene-tta-codon");
  chart.selectAll("text")
    .data(all_orfs)
  .enter().append("text")
    .attr("x", function(d) { return scale(d.start); })
    .attr("y", (single_cluster_height * i) + svgene.label_height + offset/2)
    .attr("class", "svgene-locustag")
    .attr("id", function(d) { return idx + "-cluster" + cluster.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-label"; })
    .text(function(d) { return d.locus_tag; });

};

svgene.drawUnorderedClusterOrfs = function(cluster, chart, all_orfs, scale,
                                           i, idx, height, width,xAxis,
                                           single_cluster_height, offset) {
	
  chart.append('g')
    .attr('transform', 'translate('+0+',' + 70 + ')')
    .attr('class', 'axis')
    .call(xAxis);				   
											   
  chart.selectAll("rect")
    .data(all_orfs)
  .enter().append("rect")
    .attr("x", function(d) { return scale(d.start);})
    .attr("y", (single_cluster_height * i) + svgene.label_height + offset)
    .attr("height", height - (2 * offset))
    .attr("width", function(d) { return scale(d.end) - scale(d.start)})
    .attr("opacity", 0.6)
//    .attr("rx", 3)
//    .attr("ry", 3)
    .attr("class", function(d) { return "svgene-type-" + d.type + " svgene-orf"; })
    .attr("id", function(d) { return idx + "-cluster" + cluster.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-orf"; })
    .attr("style", function(d) { if (d.color !== undefined) { return "fill:" + d.color; } })
  chart.selectAll("text")
    .data(all_orfs)
  .enter().append("text")
    .attr("x", function(d) { return scale(d.start); })
    .attr("y", (single_cluster_height * i) + svgene.label_height + offset/2)
    .attr("class", "svgene-locustag")
    .attr("id", function(d) { return idx + "-cluster" + cluster.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-label"; })
//    .text(function(d) { return d.locus_tag; });
};

svgene.drawClusters = function(id, clusters, height, width) {
  var container = d3.select("#" + id);
  var single_cluster_height = height + (2 * svgene.label_height);
  container.selectAll("svg").remove();
  container.selectAll("div").remove();
  var chart = container.append("svg")
    .attr("height", single_cluster_height * clusters.length)
    .attr("width", width + svgene.extra_label_width);
  var all_orfs = [];
  var all_borders = [];
  var all_ttas = [];

  var svg = d3.select("#"+id+"_l")

// Handmade legend

svg.append("rect").attr("x",410).attr("y",10).attr("r", 6).style("fill", "blue").attr("width",16).attr("height",16).attr("stroke","black")
svg.append("rect").attr("x",525).attr("y",10).attr("r", 6).style("fill", "lightpink").attr("width",16).attr("height",16).attr("stroke","black")
svg.append("rect").attr("x",600).attr("y",10).attr("r", 6).style("fill", "red").attr("width",16).attr("height",16).attr("stroke","black")
//svg.append("rect").attr("x",760).attr("y",10).attr("r", 6).style("fill", "lavender").attr("width",16).attr("height",16).attr("stroke","black")
svg.append("rect").attr("x",760).attr("y",10).attr("r", 6).style("fill", "#00B050").attr("width",16).attr("height",16).attr("stroke","black")
//svg.append("rect").attr("x",895).attr("y",10).attr("r", 6).style("fill", "greenyellow").attr("width",16).attr("height",16).attr("stroke","black")
svg.append("rect").attr("x",895).attr("y",10).attr("r", 6).style("fill", "#C0C0C0").attr("width",16).attr("height",16).attr("stroke","black")

svg.append("rect").attr("x",410).attr("y",35).attr("r", 6).style("fill", "brown").attr("width",16).attr("height",16).attr("stroke","black")
svg.append("rect").attr("x",525).attr("y",35).attr("r", 6).style("fill", "orange").attr("width",16).attr("height",16).attr("stroke","black")
svg.append("rect").attr("x",600).attr("y",35).attr("r", 6).style("fill", "#ba8448").attr("width",16).attr("height",16).attr("stroke","black")
svg.append("rect").attr("x",760).attr("y",35).attr("r", 6).style("fill", "#640B0F").attr("width",16).attr("height",16).attr("stroke","black")


svg.append("rect").attr("x",410).attr("y",60).attr("r", 6).style("fill", "gray").attr("width",16).attr("height",16).attr("stroke","black")
svg.append("rect").attr("x",525).attr("y",60).attr("r", 6).style("fill", "yellow").attr("width",16).attr("height",16).attr("stroke","black")
//svg.append("rect").attr("x",600).attr("y",60).attr("r", 6).style("fill", "#00B050").attr("width",16).attr("height",16).attr("stroke","black")
svg.append("rect").attr("x",600).attr("y",60).attr("r", 6).style("fill", "#03A89E").attr("width",16).attr("height",16).attr("stroke","black")
svg.append("rect").attr("x",760).attr("y",60).attr("r", 6).style("fill", "#FFFFCD").attr("width",16).attr("height",16).attr("stroke","black")
//svg.append("rect").attr("x",895).attr("y",60).attr("r", 6).style("fill", "lightgreen").attr("width",16).attr("height",16).attr("stroke","black")
svg.append("rect").attr("x",895).attr("y",60).attr("r", 6).style("fill", "#DCDCDC").attr("width",16).attr("height",16).attr("stroke","black")


svg.append("text").attr("x", 430).attr("y", 20).text("Integrase").style("font-size", "15px").attr("alignment-baseline","middle")
svg.append("text").attr("x", 545).attr("y", 20).text("T4SS").style("font-size", "15px").attr("alignment-baseline","middle")
svg.append("text").attr("x", 620).attr("y", 20).text("Antibiotic resistance").style("font-size", "15px").attr("alignment-baseline","middle")
svg.append("text").attr("x", 780).attr("y", 20).text("Defense system").style("font-size", "15px").attr("alignment-baseline","middle")
svg.append("text").attr("x", 915).attr("y", 20).text("Others").style("font-size", "15px").attr("alignment-baseline","middle")

svg.append("text").attr("x", 430).attr("y", 45).text("Relaxase").style("font-size", "15px").attr("alignment-baseline","middle")
svg.append("text").attr("x", 545).attr("y", 45).text("T4CP").style("font-size", "15px").attr("alignment-baseline","middle")
svg.append("text").attr("x", 620).attr("y", 45).text("Virulence factor").style("font-size", "15px").attr("alignment-baseline","middle")
svg.append("text").attr("x", 780).attr("y", 45).text("Degradation").style("font-size", "15px").attr("alignment-baseline","middle")


svg.append("text").attr("x", 430).attr("y", 70).text("Flank region").style("font-size", "15px").attr("alignment-baseline","middle")
svg.append("text").attr("x", 545).attr("y", 70).text("IS/Tn").style("font-size", "15px").attr("alignment-baseline","middle")
svg.append("text").attr("x", 620).attr("y", 70).text("Metal resistace").style("font-size", "15px").attr("alignment-baseline","middle")
svg.append("text").attr("x", 780).attr("y", 70).text("Symbiosis").style("font-size", "15px").attr("alignment-baseline","middle")
svg.append("text").attr("x", 915).attr("y", 70).text("Hyp").style("font-size", "15px").attr("alignment-baseline","middle")


  for (i=0; i < clusters.length; i++) {
      var cluster = clusters[i];
      all_orfs.push.apply(all_orfs, cluster.orfs.sort(svgene.sort_biosynthetic_orfs_last));
      all_borders.push.apply(all_borders, cluster.borders ? cluster.borders : []);
      all_ttas.push.apply(all_ttas, cluster.tta_codons ? cluster.tta_codons: []);
      var idx = svgene.unique_id++;
      var offset = height/10;
      var x = d3.scale.linear()
        .domain([cluster.start, cluster.end])
        .range([0, width]);
		
      var xAxis = d3.svg.axis()
     .scale(x)      //指定比例尺
     .orient("bottom")   //指定刻度的方向
     .ticks(5); 
		
      if (cluster.unordered) {
          svgene.drawUnorderedClusterOrfs(cluster, chart, all_orfs, x,
                                          i, idx, height, width,xAxis,
                                          single_cluster_height, offset);
      } else {
          svgene.drawOrderedClusterOrfs(cluster, chart, all_orfs, all_borders, all_ttas,
                                        x, i, idx, height, width,xAxis,
                                        single_cluster_height, offset);
      }
      container.selectAll("div")
        .data(all_orfs)
        .enter().append("div")
        .attr("class", "svgene-tooltip")
        .attr("id", function(d) { return idx + "-cluster" + cluster.idx + "-" + svgene.tag_to_id(d.locus_tag) + "-tooltip"; })
        .html(function(d) { return d.description});
  }
  for (i=0; i < clusters.length; i++) {
      var cluster = clusters[i];
      if (cluster.label !== undefined) {
        chart.append("text")
            .text(cluster.label)
            .attr("class", "svgene-clusterlabel")
            .attr("x", function() { return width + svgene.extra_label_width - this.getComputedTextLength() - 5})
            .attr("y", function() { return (single_cluster_height * i) + svgene.label_height } )
            .attr("font-size", svgene.label_height);
      }
  }
  svgene.init();
};

svgene.sort_biosynthetic_orfs_last = function(a, b) {
    if ((a.type != "biosynthetic" && b.type != "biosynthetic") ||
       (a.type == "biosynthetic" && b.type == "biosynthetic")) {
        return a.start - b.start;
    };
    if (a.type == "biosynthetic") {
       return 1;
    }
    return -1;
};

svgene.tag_to_id = function(tag) {
    return tag.replace(/(:|\.)/g, '-').replace(/-orf/g, '_orf');
}


svgene.tooltip_handler = function(ev) {
    var id = $(this).attr("id").replace("-orf", "-tooltip");
    var tooltip = $("#"+id);

    if (svgene.active_tooltip) {
        svgene.active_tooltip.hide();
    }
    svgene.active_tooltip = tooltip;

    if (tooltip.css("display") == 'none') {
        var offset = $(this).offset();
        tooltip.css("middle", offset.left + 10);
        var this_parent = $(this).parent();
        var top_offset = this_parent.height()/(this_parent.children('line').length * 2);
        tooltip.css("bottom", offset.top + top_offset);
        tooltip.show();
        tooltip.click(function(){$(this).hide()});
        var timeout = setTimeout(function(){ tooltip.slideUp("fast") }, 5000);
        tooltip.data("timeout", timeout);
        tooltip.mouseover(function() {
            clearTimeout(tooltip.data("timeout"));
        }).mouseout(function() {
            timeout = setTimeout(function(){ tooltip.slideUp("fast") }, 5000);
            tooltip.data("timeout", timeout);
        });
    } else {
        tooltip.hide();
    }
};

svgene.init = function() {
    $(".svgene-orf").mouseover(function(e) {
        var id = $(this).attr("id").replace("-orf", "-label");
		this.style.opacity = 1;
        $("#"+id).show();
    }).mouseout(function(e) {
        var id = $(this).attr("id").replace("-orf", "-label");
		this.style.opacity = 0.8;
        $("#"+id).hide();
    }).click(svgene.tooltip_handler);
    $(".svgene-textarea").click(function(event) {
        event.stopPropagation();
    });
};
