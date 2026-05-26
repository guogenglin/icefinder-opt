
let container = document.getElementById('gcmap'),
  highlightPoint = function(e) {
    var chart,
      point,
      i,
      event;
    for (i = 0; i < Highcharts.charts.length; i = i + 1) {
      chart = Highcharts.charts[i];
      // Find coordinates within the chart
      event = chart.pointer.normalize(e);
      // Get the hovered point
      point = chart.series[0].searchPoint(event, true);
      if (point) {
        point.highlight(e);
      }
    }
  };
['mousemove', 'touchmove', 'touchstart'].forEach(function (eventType) {
  container.addEventListener(
    eventType,
    function (e) {
      highlightPoint(e);
    }
  );
});
/**
 * 重写内部的方法， 这里是将提示框即十字准星的隐藏函数关闭
 */
Highcharts.Pointer.prototype.reset = function () {
  return undefined;
};
/**
 * 高亮当前的数据点，并设置鼠标滑动状态及绘制十字准星线
 */
Highcharts.Point.prototype.highlight = function (event) {
  this.onMouseOver(); // 显示鼠标激活标识
  this.series.chart.tooltip.refresh(this); // 显示提示框
  this.series.chart.xAxis[0].drawCrosshair(event, this); // 显示十字准星线
};
/**
 * 同步缩放效果，即当一个图表进行了缩放效果，其他图表也进行缩放
 */
function syncExtremes(e) {
  var thisChart = this.chart;
  if (e.trigger !== 'syncExtremes') { 
    Highcharts.each(Highcharts.charts, function (chart) {
      if (chart !== thisChart) {
        if (chart.xAxis[0].setExtremes) { 
          chart.xAxis[0].setExtremes(e.min, e.max, undefined, false, { trigger: 'syncExtremes' });
        }
      }
    });
  }
}
// 使用定义的JSON数据
var activity = jsonData;

activity.datasets.forEach(function (dataset, i) {
  // Add X values
  dataset.data = Highcharts.map(dataset.data, function (val, j) {
    return [activity.xData[j], val];
  });
  var chartDiv = document.createElement('div');
  chartDiv.className = 'chart';
  container.appendChild(chartDiv);
  Highcharts.chart(chartDiv, {
    chart: {
      width: 990,
      height: 200,
      marginLeft: 45, // Keep all charts left aligned
      spacingTop: 20,
      spacingBottom: 20
    },
    title: {
      text: dataset.name,
      align: 'left',
      margin: 0,
      x: 30
    },
    credits: {
      enabled: false
    },
    legend: {
      enabled: false
    },
    xAxis: {
      crosshair: true,
      events: {
        setExtremes: syncExtremes
      },
      labels: {
        format: '{value} k'
      }
    },
    yAxis: {
      title: {
        text: null
      }
    },
    tooltip: {
      positioner: function () {
        return {
          // right aligned
          x: this.chart.chartWidth - this.label.width,
          y: 10 // align to title
        };
      },
      borderWidth: 0,
      backgroundColor: 'none',
      pointFormat: '{point.y}',
      headerFormat: '',
      shadow: false,
      style: {
        fontSize: '12px'
      },
      valueDecimals: dataset.valueDecimals
    },
    series: [{
      data: dataset.data,
      name: dataset.name,
      type: dataset.type,
      lineWidth:1,
      color: Highcharts.getOptions().colors[i],
      fillOpacity: 0.3,
      tooltip: {
        valueSuffix: ' ' + dataset.unit
      }
    }]
  });
});
