#!/usr/bin/miniconda3/bin/python3
"""
---------------------------------------------
Generate HTML report from NIPT JSON file

Author: Hyukjung Kwon
Updated: 2025-06-05
---
"""

import os
import sys
import json
import argparse
import logging
from pathlib import Path
from datetime import datetime

__author__ = 'Kenneth Kwon'
__email__ = "joykwon77@gmail.com"
__version__ = '1.0'

logger = logging.getLogger(__name__)

def file_check(parser, arg):
    if not os.path.exists(arg):
        parser.error(f"The file {arg} does not exist!")
    else:
        return str(arg)

def generate_html_template():
    """Generate HTML template with CSS styling"""
    return """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NIPT Analysis Report - {sample_id}</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background-color: #f5f5f5;
            padding: 20px;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
            position: relative;
        }}
        
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        
        .header .subtitle {{
            font-size: 1.2em;
            opacity: 0.9;
        }}
        
        .algorithm-version {{
            position: absolute;
            top: 20px;
            right: 30px;
            background: rgba(255,255,255,0.2);
            padding: 8px 15px;
            border-radius: 20px;
            font-size: 0.9em;
        }}
        
        .content {{
            padding: 30px;
        }}
        
        .section {{
            margin-bottom: 40px;
            background: #fff;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }}
        
        .section-header {{
            background: #f8f9fa;
            padding: 20px;
            border-bottom: 3px solid #667eea;
            border-radius: 8px 8px 0 0;
        }}
        
        .section-title {{
            font-size: 1.8em;
            color: #2c3e50;
            margin: 0;
            display: flex;
            align-items: center;
        }}
        
        .section-title::before {{
            content: "";
            width: 4px;
            height: 30px;
            background: #667eea;
            margin-right: 15px;
            border-radius: 2px;
        }}
        
        .section-content {{
            padding: 25px;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        
        th {{
            background: #667eea;
            color: white;
            font-weight: 600;
            text-transform: uppercase;
            font-size: 0.9em;
            letter-spacing: 0.5px;
        }}
        
        tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}
        
        tr:hover {{
            background-color: #e3f2fd;
            transition: background-color 0.3s ease;
        }}
        
        .btn {{
            display: inline-block;
            padding: 8px 16px;
            margin: 2px;
            background: #667eea;
            color: white;
            text-decoration: none;
            border-radius: 5px;
            font-size: 0.9em;
            font-weight: 500;
            transition: all 0.3s ease;
            border: none;
            cursor: pointer;
        }}
        
        .btn:hover {{
            background: #5a6fd8;
            transform: translateY(-2px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.2);
        }}
        
        .btn-plot {{
            background: #28a745;
        }}
        
        .btn-plot:hover {{
            background: #218838;
        }}
        
        .btn-report {{
            background: #fd7e14;
        }}
        
        .btn-report:hover {{
            background: #e96500;
        }}
        
        .btn-image {{
            background: #6f42c1;
        }}
        
        .btn-image:hover {{
            background: #5a32a3;
        }}
        
        .status-pass {{
            color: #28a745;
            font-weight: bold;
        }}
        
        .status-fail {{
            color: #dc3545;
            font-weight: bold;
        }}
        
        .status-normal {{
            color: #28a745;
            font-weight: bold;
        }}
        
        .status-abnormal {{
            color: #dc3545;
            font-weight: bold;
        }}
        
        .high-risk {{
            color: #dc3545;
            font-weight: bold;
        }}
        
        .low-risk {{
            color: #28a745;
            font-weight: bold;
        }}
        
        .detected {{
            color: #dc3545;
            font-weight: bold;
        }}
        
        .not-detected {{
            color: #28a745;
            font-weight: bold;
        }}
        
        .info-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        
        .info-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        
        .info-label {{
            font-weight: bold;
            color: #495057;
            margin-bottom: 5px;
        }}
        
        .info-value {{
            font-size: 1.1em;
            color: #2c3e50;
        }}
        
        .button-group {{
            display: flex;
            flex-wrap: wrap;
            gap: 5px;
            margin-top: 10px;
        }}
        
        .sub-section {{
            margin: 30px 0;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 6px;
            border-left: 4px solid #667eea;
        }}
        
        .sub-title {{
            font-size: 1.3em;
            color: #2c3e50;
            margin-bottom: 15px;
            font-weight: 600;
        }}
        
        .footer {{
            text-align: center;
            padding: 20px;
            background: #f8f9fa;
            color: #666;
            font-size: 0.9em;
        }}
        
        @media (max-width: 768px) {{
            .container {{
                margin: 10px;
                border-radius: 0;
            }}
            
            .header {{
                padding: 20px;
            }}
            
            .header h1 {{
                font-size: 2em;
            }}
            
            .content {{
                padding: 15px;
            }}
            
            table {{
                font-size: 0.9em;
            }}
            
            th, td {{
                padding: 8px 10px;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="algorithm-version">Algorithm {algorithm_version}</div>
            <h1>NIPT Analysis Report</h1>
            <div class="subtitle">Order ID: {sample_id}</div>
        </div>
        
        <div class="content">
            {content}
        </div>
        
        <div class="footer">
            Generated on {timestamp} | NIPT Analysis System v{version}
        </div>
    </div>
</body>
</html>"""

def generate_file_button(file_path, button_text, button_class="btn", sample_id=None):
    """Generate file button HTML with proper path resolution"""
    if isinstance(file_path, list):
        return ""

    if file_path and file_path != "":
        # If sample_id is provided, prepend it to create relative path
        if sample_id and not file_path.startswith('http'):
            # Remove leading slash if present and prepend sample_id
            clean_path = file_path.lstrip('/')
            adjusted_path = f"{sample_id}/{clean_path}"
        else:
            adjusted_path = file_path
        return f'<a href="{adjusted_path}" target="_blank" class="{button_class}">{button_text}</a>'
    return ""

def generate_file_button_list(file_paths, button_text_prefix, button_class="btn", sample_id=None):
    """Generate multiple file buttons if file_paths is a list"""
    if isinstance(file_paths, list):
        buttons = []
        for idx, path in enumerate(file_paths):
            label = f"{button_text_prefix} {idx+1}" if len(file_paths) > 1 else button_text_prefix
            buttons.append(generate_file_button(path, label, button_class, sample_id))
        return " ".join(buttons)
    else:
        return generate_file_button(file_paths, button_text_prefix, button_class, sample_id)


def generate_lab_test_section(data):
    """Generate Lab Test section HTML"""
    lab_test = data.get('lab_test', {})
    
    html = """
    <div class="section">
        <div class="section-header">
            <h2 class="section-title">🔬 Laboratory Quality Control</h2>
        </div>
        <div class="section-content">
            <table>
                <thead>
                    <tr>
                        <th>Test Item</th>
                        <th>Result</th>
                        <th>Status</th>
                    </tr>
                </thead>
                <tbody>
    """
    
    test_items = {
        'sample_suitability': 'Sample Suitability',
        'dna_quality': 'DNA Quality',
        'library_quality': 'Library Quality', 
        'ngs_data_quality': 'NGS Data Quality',
        'reference_material_test': 'Reference Material Test'
    }
    
    for key, label in test_items.items():
        value = lab_test.get(key, 'N/A')
        status_class = 'status-pass' if value == 'Pass' else 'status-fail'
        html += f"""
                    <tr>
                        <td>{label}</td>
                        <td><span class="{status_class}">{value}</span></td>
                        <td><span class="{status_class}">{'✓' if value == 'Pass' else '✗'}</span></td>
                    </tr>
        """
    
    html += """
                </tbody>
            </table>
        </div>
    </div>
    """
    return html


def generate_final_results_section(data):
    """Generate Final Results section HTML with proper list handling"""
    final_results = data.get('final_results', {})
    
    html = """
    <div class="section">
        <div class="section-header">
            <h2 class="section-title">📊 Final Results Summary</h2>
        </div>
        <div class="section-content">
            <div class="info-grid">
    """
    
    def format_list_value(value, label):
        """리스트 값을 적절하게 포맷팅"""
        if isinstance(value, list):
            if not value:  # 빈 리스트
                if label in ('Trisomy Result', 'Microdeletion Result'):
                    return "Low Risk"
                return "Not Detected"
            elif len(value) == 1:  # 단일 항목
                return str(value[0])
            else:  # 여러 항목
                # 모든 리스트를 쉼표로 구분해서 표시
                return ', '.join(str(item) for item in value)
        else:
            return str(value) if value is not None else 'N/A'
    
    def get_result_status_class(value, label):
        """결과 값에 따른 CSS 클래스 결정"""
        # 1) 리스트일 때
        if isinstance(value, list):
            normal_terms = {'Low Risk', 'Normal', 'Not Detected'}
            # 1) 빈 리스트: 정상
            if not value:
                return "status-normal"
            # 2) 리스트에 하나라도 정상 용어가 아닌 항목이 있으면 비정상
            if any(item not in normal_terms for item in value):
                return "status-abnormal"
            # 3) 모두 정상 용어만 있을 경우: 정상
            return "status-normal"

        # 2) 스칼라 값일 때 (기존 로직)
        if 'Result' in label:
            if value in ['Low Risk', 'Normal', 'Not Detected']:
                return "status-normal"
            if value in ['High Risk', 'Abnormal', 'Detected']:
                return "status-abnormal"

        if label.lower() == 'sample_bias' and value == 'PASS':
            return "status-pass"

        return "" 

    # 정보 항목들 - 정확한 JSON 키 이름 사용
    info_items = [
        ('Order ID', final_results.get('order_id', 'N/A')),
        ('Fetal Fraction (YFF)', final_results.get('fetal_fraction_yff', 'N/A')),
        ('Fetal Fraction (seqFF)', final_results.get('fetal_fraction_seqff', 'N/A')),
        ('FF Ratio', final_results.get('ff_ratio', 'N/A')),
        ('Sample Bias', final_results.get('sample_bias', 'N/A')),
        ('Fetal Gender', final_results.get('fetal_gender', 'N/A')),
        ('Trisomy Result', final_results.get('trisomy_result', [])),  # final_ 없음!
        ('Microdeletion Result', final_results.get('md_result', []))  # final_ 없음!
    ]
    
    for label, value in info_items:
        # 값 포맷팅
        formatted_value = format_list_value(value, label)
        
        # CSS 클래스 결정
        value_class = get_result_status_class(value, label)
        
        html += f"""
                <div class="info-card">
                    <div class="info-label">{label}</div>
                    <div class="info-value"><span class="{value_class}">{formatted_value}</span></div>
                </div>
        """
    
    html += """
            </div>
        </div>
    </div>
    """
    
    return html

# 사용 예시
if __name__ == "__main__":
    # 테스트 데이터
    test_data = {
        'final_results': {
            'order_id': '2506040012',
            'fetal_fraction_yff': 'N/A',
            'fetal_fraction_seqff': '2.44',
            'ff_ratio': '0.28',
            'sample_bias': 'PASS',
            'fetal_gender': 'Female',
            'final_trisomy_result': ['Trisomy21'],
            'final_md_result': ['md8', 'other_md320', 'other_md87']
        }
    }
    
    logger.info("=== 옵션 1: 배지 스타일 (추천) ===")
    result1 = generate_final_results_section(test_data)
    logger.info("HTML 생성 완료")


def generate_trisomy_results_section(data):
    """Generate Trisomy Results section HTML"""
    trisomy_results = data.get('trisomy_results', [])
    
    html = """
    <div class="section">
        <div class="section-header">
            <h2 class="section-title">🧬 Trisomy Analysis Results</h2>
        </div>
        <div class="section-content">
            <table>
                <thead>
                    <tr>
                        <th>Item</th>
                        <th>Disease Name</th>
                        <th>Result</th>
                        <th>Risk Before (Single)</th>
                        <th>Risk Before (Twin)</th>
                        <th>Risk After</th>
                        <th>PPV</th>
                        <th>NPV</th>
                    </tr>
                </thead>
                <tbody>
    """
    
    for result in trisomy_results:
        result_value = result.get('result', 'N/A')
        result_class = "low-risk" if result_value == 'Low Risk' else "high-risk"
        
        html += f"""
                    <tr>
                        <td><strong>{result.get('item', 'N/A')}</strong></td>
                        <td>{result.get('disease_name', 'N/A')}</td>
                        <td><span class="{result_class}">{result_value}</span></td>
                        <td>{result.get('risk_before_single') or 'N/A'}</td>
                        <td>{result.get('risk_before_twin') or 'N/A'}</td>
                        <td>{result.get('risk_after') or 'N/A'}</td>
                        <td>{result.get('ppv') or 'N/A'}</td>
                        <td>{result.get('npv') or 'N/A'}</td>
                    </tr>
        """
    
    html += """
                </tbody>
            </table>
        </div>
    </div>
    """
    return html

def generate_trisomy_details_section(data, sample_id=None):
    """Generate Trisomy Details section HTML"""
    trisomy_details = data.get('trisomy_details', {})
    
    html = """
    <div class="section">
        <div class="section-header">
            <h2 class="section-title">📈 Chromosome Analysis Details</h2>
        </div>
        <div class="section-content">
    """
    
    data_sources = [
        ('orig', 'Original Analysis'),
        ('fetus', 'Fetus Analysis'),
        ('mom', 'Maternal Analysis')
    ]
    
    for data_src, title in data_sources:
        if data_src not in trisomy_details:
            continue
            
        details = trisomy_details[data_src]
        result_table = details.get('result_table', {})
        
        html += f"""
            <div class="sub-section">
                <h3 class="sub-title">{title}</h3>
                
                <div class="button-group">
        """
        
        # Add file buttons
        file_buttons = [
            (f'{data_src}_ezd_plot', "EZD Plot", "btn btn-plot"),
            (f'{data_src}_prizm_chr_plot', "PRIZM Chr Plot", "btn btn-plot"),
            (f'{data_src}_prizm_10mb_plot', "PRIZM 10MB Plot", "btn btn-plot"),
            (f'{data_src}_wc_plot', "WC Plot", "btn btn-plot"),
            (f'{data_src}_wc_result', "WC Report", "btn btn-report"),
            (f'{data_src}_wcx_plot', "WCX Plot", "btn btn-image"),
            (f'{data_src}_wcx_result', "WCX Report", "btn btn-report")
        ]
        
        for file_key, button_text, button_class in file_buttons:
            file_path = details.get(file_key)
            if file_path:
                html += generate_file_button(file_path, button_text, button_class, sample_id)
        
        html += """
                </div>
                
                <table>
                    <thead>
                        <tr>
                            <th>Chromosome</th>
                            <th>EZD Detection</th>
                            <th>PRIZM Detection</th>
                            <th>Z-Score</th>
                            <th>UAR (%)</th>
                            <th>Z-Score Threshold</th>
                            <th>UAR Threshold</th>
                        </tr>
                    </thead>
                    <tbody>
        """
        def chromosome_sort_key(chr_item):
            """염색체를 올바른 순서로 정렬하기 위한 키 함수"""
            chr_key = chr_item[0]
            
            # 'Chromosome' 제거하고 실제 염색체 이름만 추출
            if chr_key.startswith('Chromosome '):
                chr_name = chr_key.replace('Chromosome ', '').strip()
            else:
                chr_name = chr_key
            
            # 숫자 염색체는 정수로 변환하여 정렬
            if chr_name.isdigit():
                return (0, int(chr_name))  # (타입, 번호)
            # X 염색체
            elif chr_name == 'X':
                return (1, 0)  # 숫자 염색체 다음
            # Y 염색체  
            elif chr_name == 'Y':
                return (1, 1)  # X 염색체 다음
            # 기타
            else:
                return (2, chr_name)  # 맨 마지막

        # Sort chromosomes for better display
        sorted_chromosomes = sorted(result_table.items(), key=chromosome_sort_key)
        
        for chr_key, chr_data in sorted_chromosomes:
            ezd_detection = chr_data.get('EZD Detection', 'N/A')
            ezd_detection_class = "low-risk" if ezd_detection == 'Low Risk' else ("high-risk" if ezd_detection == 'High Risk' else "low-risk")
            
            prizm_detection = chr_data.get('PRIZM Detection', 'N/A')
            prizm_detection_class = "low-risk" if prizm_detection == 'Low Risk' else ("high-risk" if prizm_detection == 'High Risk' else "low-risk")

            html += f"""
                        <tr>
                            <td><strong>{chr_key}</strong></td>
                            <td><span class="{ezd_detection_class}">{ezd_detection}</span></td>
                            <td><span class="{prizm_detection_class}">{prizm_detection}</span></td>
                            <td>{chr_data.get('Z-score', 'N/A')}</td>
                            <td>{chr_data.get('UAR(%)', 'N/A')}</td>
                            <td>{chr_data.get('Z-score threshold', 'N/A')}</td>
                            <td>{chr_data.get('UAR threshold', 'N/A')}</td>
                        </tr>
            """
        
        html += """
                    </tbody>
                </table>
            </div>
        """
    
    html += """
        </div>
    </div>
    """
    return html

def generate_microdeletion_section(data, sample_id=None):
    """Generate Microdeletion section HTML"""
    md_results = data.get('md_results', {})
    md_details = data.get('md_details', {})
    
    html = """
    <div class="section">
        <div class="section-header">
            <h2 class="section-title">🔍 Microdeletion Analysis</h2>
        </div>
        <div class="section-content">
    """
    
    # MD Results Summary
    html += """
            <div class="sub-section">
                <h3 class="sub-title">Summary Results</h3>
                <table>
                    <thead>
                        <tr>
                            <th>Item</th>
                            <th>Location</th>
                            <th>Disease Name</th>
                            <th>Result</th>
                        </tr>
                    </thead>
                    <tbody>
    """
    
    result_table = md_results.get('result_table', [])
    for result in result_table:
        result_value = result.get('result', 'N/A')
        result_class = "high-risk" if result_value == 'High Risk' else "low-risk"
        html += f"""
                        <tr>
                            <td><strong>{result.get('item', 'N/A')}</strong></td>
                            <td>{result.get('location', 'N/A')}</td>
                            <td>{result.get('disease_name', 'N/A')}</td>
                            <td><span class="{result_class}">{result_value}</span></td>
                        </tr>
        """
    
    html += """
                    </tbody>
                </table>
            </div>
    """
    
    # MD Details
    md_sections = [
        ('md8_results', '8 Common Syndromes'),
        ('md108_results', '108 Extended Syndromes'),
        ('md320_results', '320 Extended Syndromes'),
        ('md87_results', '87 Extended Syndromes'),
        ('md141_results', '141 Extended Syndromes')
    ]
    
    for md_section, section_title in md_sections:
        if md_section not in md_details:
            continue
            
        section_data = md_details[md_section]
        
        html += f"""
            <div class="sub-section">
                <h3 class="sub-title">{section_title}</h3>
        """
        
        data_sources = [
            ('orig', 'Original'),
            ('fetus', 'Fetus'),
            ('mom', 'Maternal')
        ]
        
        for data_src, src_title in data_sources:
            if data_src not in section_data:
                continue
                
            src_data = section_data[data_src]
            
            html += f"""
                <h4>{src_title} Analysis</h4>
                <div class="button-group">
            """
            
            # Image buttons
            image_data = src_data.get('image', {})
            if image_data and isinstance(image_data, dict):
                html += generate_file_button(image_data.get('WC'), 'WC Plot', 'btn btn-plot', sample_id)
                html += generate_file_button(image_data.get('WCX'), 'WCX Plot', 'btn btn-image', sample_id)
            
            html += """
                </div>
            """
            
            # MD items with detection
            detected_items = []
            for key, value in src_data.items():
                if key.startswith('md') and isinstance(value, dict):
                    disease_name = value.get('disease_name')
                    if disease_name and disease_name not in [None, 'N/A', '-', '']:
                        detected_items.append((key, value))
            
            if detected_items:
                html += """
                <table>
                    <thead>
                        <tr>
                            <th>MD Item</th>
                            <th>Disease Name</th>
                            <th>Target Region</th>
                            <th>WC Detection</th>
                            <th>WCX Detection</th>
                            <th>WC Details</th>
                            <th>WCX Details</th>
                        </tr>
                    </thead>
                    <tbody>
                """
                
                for key, value in detected_items:
                    disease_name = value.get('disease_name', 'N/A')
                    target_region = value.get('target_region', 'N/A')
                    
                    detection = value.get('detection', {})
                    wc_detection = detection.get('WC') or 'Low Risk'
                    wcx_detection = detection.get('WCX') or 'Low Risk'
                    
                    wc_class = "high-risk" if wc_detection == 'High Risk' else "low-risk"
                    wcx_class = "high-risk" if wcx_detection == 'High Risk' else "low-risk"
                    
                    # WC Details
                    wc_details = ""
                    if wc_detection == 'High Risk':
                        detected_region = value.get('detected_region', {}).get('WC', '')
                        z_score = value.get('z_score', {}).get('WC', '')
                        length = value.get('length', {}).get('WC', '')
                        link = value.get('detected_region_link', {}).get('WC', '')
                        
                        wc_details = f"Region: {detected_region}<br>Z-score: {z_score}<br>Length: {length}"
                        if link:
                            wc_details += f'<br><a href="{link}" target="_blank" class="btn" style="font-size:0.8em;">Decipher</a>'
                    
                    # WCX Details
                    wcx_details = ""
                    if wcx_detection == 'High Risk':
                        detected_region = value.get('detected_region', {}).get('WCX', '')
                        z_score = value.get('z_score', {}).get('WCX', '')
                        length = value.get('length', {}).get('WCX', '')
                        link = value.get('detected_region_link', {}).get('WCX', '')
                        image = value.get('image', {}).get('WCX', '')
                        
                        wcx_details = f"Region: {detected_region}<br>Z-score: {z_score}<br>Length: {length}"
                        if link:
                            wcx_details += f'<br><a href="{link}" target="_blank" class="btn" style="font-size:0.8em;">Decipher</a>'
                        if image:
                            # Kenneth : 250719 : this part comes as list type. startswith makes error
                            image_button = generate_file_button_list(image, "Image", "btn btn-image", sample_id)
                            if image_button:
                                wcx_details += f'<br><div style="font-size:0.8em;">{image_button}</div>'

                            #wcx_details += f'<br><a href="{image}" target="_blank" class="btn btn-image" style="font-size:0.8em;">Image</a>'
                    
                    html += f"""
                        <tr>
                            <td><strong>{key}</strong></td>
                            <td>{disease_name}</td>
                            <td>{target_region}</td>
                            <td><span class="{wc_class}">{wc_detection}</span></td>
                            <td><span class="{wcx_class}">{wcx_detection}</span></td>
                            <td>{wc_details}</td>
                            <td>{wcx_details}</td>
                        </tr>
                    """
                
                html += """
                    </tbody>
                </table>
                """
        
        html += """
            </div>
        """
    
    html += """
        </div>
    </div>
    """
    return html

def generate_quality_control_section(data, sample_id=None):
    """Generate Quality Control section HTML"""
    # 올바른 데이터 접근: NIPT 안의 quality_control
    nipt_data = data.get('NIPT', {})
    qc_data = nipt_data.get('quality_control', {})

    html = """
    <div class="section">
        <div class="section-header">
            <h2 class="section-title">🎯 Quality Control Results</h2>
        </div>
        <div class="section-content">
    """

    # Sequencing Metrics
    sequencing_metrics = qc_data.get('sequencing_metrics', {})
    if sequencing_metrics:
        html += """
            <div class="sub-section">
                <h3 class="sub-title">Sequencing Metrics</h3>
                <table>
                    <thead>
                        <tr>
                            <th>Metric</th>
                            <th>Value</th>
                            <th>Unit</th>
                            <th>Threshold</th>
                            <th>Status</th>
                        </tr>
                    </thead>
                    <tbody>
        """

        metric_labels = {
            'total_reads': 'Total Reads',
            'mapped_reads': 'Mapped Reads',
            'mapping_rate': 'Mapping Rate',
            'duplication_rate': 'Duplication Rate',
            'mean_mapping_quality': 'Mean Mapping Quality',
            'mean_coverage': 'Mean Coverage',
            'gc_content': 'GC Content'
        }

        for key, label in metric_labels.items():
            if key in sequencing_metrics:
                metric = sequencing_metrics[key]
                value = metric.get('value', 'N/A')
                unit = metric.get('unit', '')
                threshold = metric.get('threshold', 'N/A')
                status = metric.get('status', 'N/A')
                status_class = 'status-pass' if status == 'PASS' else 'status-fail'

                # Format large numbers - 수정된 부분
                try:
                    # 문자열인 경우 (예: "0.3662X") 그대로 사용
                    if isinstance(value, str):
                        formatted_value = value
                    elif isinstance(value, (int, float)) and value > 1000:
                        if value > 1000000:
                            formatted_value = f"{value/1000000:.1f}M"
                        else:
                            formatted_value = f"{value/1000:.1f}K"
                    else:
                        formatted_value = str(value)
                except Exception:
                    formatted_value = str(value)

                html += f"""
                        <tr>
                            <td>{label}</td>
                            <td>{formatted_value}</td>
                            <td>{unit}</td>
                            <td>{threshold}</td>
                            <td><span class="{status_class}">{status}</span></td>
                        </tr>
                """

        html += """
                    </tbody>
                </table>
            </div>
        """

    # Analysis QC
    analysis_qc = qc_data.get('analysis_qc', {})
    if analysis_qc:
        html += """
            <div class="sub-section">
                <h3 class="sub-title">Analysis Quality Control</h3>
                <table>
                    <thead>
                        <tr>
                            <th>Parameter</th>
                            <th>Value</th>
                            <th>Unit</th>
                            <th>Threshold</th>
                            <th>Status</th>
                        </tr>
                    </thead>
                    <tbody>
        """

        analysis_labels = {
            'fetal_fraction_yff': 'Fetal Fraction (YFF)',
            'fetal_fraction_seqff': 'Fetal Fraction (seqFF)',
            'ff_ratio': 'FF Ratio',
            'sample_bias_qc': 'Sample Bias QC'
        }

        for key, label in analysis_labels.items():
            if key in analysis_qc:
                param = analysis_qc[key]
                value = param.get('value', 'N/A')
                unit = param.get('unit', '')
                threshold = param.get('threshold', 'N/A')
                status = param.get('status', 'N/A')
                status_class = 'status-pass' if status == 'PASS' else 'status-fail'

                html += f"""
                        <tr>
                            <td>{label}</td>
                            <td>{value}</td>
                            <td>{unit}</td>
                            <td>{threshold}</td>
                            <td><span class="{status_class}">{status}</span></td>
                        </tr>
                """

        html += """
                    </tbody>
                </table>
            </div>
        """

    # QC Files
    qc_files = qc_data.get('qc_files', {})
    if qc_files:
        html += """
            <div class="sub-section">
                <h3 class="sub-title">QC Reports</h3>
                <div class="button-group">
        """

        file_buttons = [
            ('Fastqc_R1_report', 'FastQC R1 Report', 'btn btn-report'),
            ('Fastqc_R2_report', 'FastQC R2 Report', 'btn btn-report'),
            ('Qualimap_report', 'Qualimap Report', 'btn btn-report')
        ]

        for file_key, button_text, button_class in file_buttons:
            file_path = qc_files.get(file_key)
            if file_path:
                html += generate_file_button(file_path, button_text, button_class, sample_id)

        html += """
                </div>
            </div>
        """

    html += """
        </div>
    </div>
    """
    return html

def generate_quality_control_section_old(data, sample_id=None):
    """Generate Quality Control section HTML"""
    qc_data = data.get('quality_control', {})
    
    html = """
    <div class="section">
        <div class="section-header">
            <h2 class="section-title">🎯 Quality Control Results</h2>
        </div>
        <div class="section-content">
    """
    
    # Sequencing Metrics
    sequencing_metrics = qc_data.get('sequencing_metrics', {})
    if sequencing_metrics:
        html += """
            <div class="sub-section">
                <h3 class="sub-title">Sequencing Metrics</h3>
                <table>
                    <thead>
                        <tr>
                            <th>Metric</th>
                            <th>Value</th>
                            <th>Unit</th>
                            <th>Threshold</th>
                            <th>Status</th>
                        </tr>
                    </thead>
                    <tbody>
        """
        
        metric_labels = {
            'total_reads': 'Total Reads',
            'mapped_reads': 'Mapped Reads',
            'mapping_rate': 'Mapping Rate',
            'duplication_rate': 'Duplication Rate',
            'mean_mapping_quality': 'Mean Mapping Quality',
            'mean_coverage': 'Mean Coverage',
            'gc_content': 'GC Content'
        }
        
        for key, label in metric_labels.items():
            if key in sequencing_metrics:
                metric = sequencing_metrics[key]
                value = metric.get('value', 'N/A')
                unit = metric.get('unit', '')
                threshold = metric.get('threshold', 'N/A')
                status = metric.get('status', 'N/A')
                status_class = 'status-pass' if status == 'PASS' else 'status-fail'
                
                # Format large numbers
                if isinstance(value, (int, float)) and value > 1000:
                    if value > 1000000:
                        formatted_value = f"{value/1000000:.1f}M"
                    else:
                        formatted_value = f"{value/1000:.1f}K"
                else:
                    formatted_value = str(value)
                
                html += f"""
                        <tr>
                            <td>{label}</td>
                            <td>{formatted_value}</td>
                            <td>{unit}</td>
                            <td>{threshold}</td>
                            <td><span class="{status_class}">{status}</span></td>
                        </tr>
                """
        
        html += """
                    </tbody>
                </table>
            </div>
        """
    
    # Analysis QC
    analysis_qc = qc_data.get('analysis_qc', {})
    if analysis_qc:
        html += """
            <div class="sub-section">
                <h3 class="sub-title">Analysis Quality Control</h3>
                <table>
                    <thead>
                        <tr>
                            <th>Parameter</th>
                            <th>Value</th>
                            <th>Unit</th>
                            <th>Threshold</th>
                            <th>Status</th>
                        </tr>
                    </thead>
                    <tbody>
        """
        
        analysis_labels = {
            'fetal_fraction_yff': 'Fetal Fraction (YFF)',
            'fetal_fraction_seqff': 'Fetal Fraction (seqFF)',
            'ff_ratio': 'FF Ratio',
            'sample_bias_qc': 'Sample Bias QC'
        }
        
        for key, label in analysis_labels.items():
            if key in analysis_qc:
                param = analysis_qc[key]
                value = param.get('value', 'N/A')
                unit = param.get('unit', '')
                threshold = param.get('threshold', 'N/A')
                status = param.get('status', 'N/A')
                status_class = 'status-pass' if status == 'PASS' else 'status-fail'
                
                html += f"""
                        <tr>
                            <td>{label}</td>
                            <td>{value}</td>
                            <td>{unit}</td>
                            <td>{threshold}</td>
                            <td><span class="{status_class}">{status}</span></td>
                        </tr>
                """
        
        html += """
                    </tbody>
                </table>
            </div>
        """
    
    # QC Files
    qc_files = qc_data.get('qc_files', {})
    if qc_files:
        html += """
            <div class="sub-section">
                <h3 class="sub-title">QC Reports</h3>
                <div class="button-group">
        """
        
        file_buttons = [
            ('Fastqc_R1_report', 'FastQC R1 Report', 'btn btn-report'),
            ('Fastqc_R2_report', 'FastQC R2 Report', 'btn btn-report'),
            ('Qualimap_report', 'Qualimap Report', 'btn btn-report')
        ]
        
        for file_key, button_text, button_class in file_buttons:
            file_path = qc_files.get(file_key)
            if file_path:
                html += generate_file_button(file_path, button_text, button_class, sample_id)
        
        html += """
                </div>
            </div>
        """
    
    html += """
        </div>
    </div>
    """
    return html

def generate_html_report(json_file, output_dir):
    """Generate complete HTML report from JSON file"""
    
    # Read JSON data
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    nipt_data = data.get('NIPT', {})
    sample_id = nipt_data.get('final_results', {}).get('order_id', 'Unknown')
    algorithm_version = nipt_data.get('algorithm_version', 'v1.0')
    
    # Check if sample directory exists and inform user
    #sample_dir = os.path.join(os.path.dirname(json_file), sample_id)
    sample_dir = os.path.dirname(json_file)
    if not os.path.exists(sample_dir):
        logger.info(f"⚠️  Warning: Sample directory '{sample_dir}' not found.")
        logger.info(f"📂 Expected structure: {sample_id}/Output_EZD/, {sample_id}/Output_PRIZM/, etc.")
        logger.info("🔗 File links in HTML may not work unless the directory structure is correct.")
    
    # Generate content sections
    content_sections = []
    
    # 1. Lab Test Section
    content_sections.append(generate_lab_test_section(nipt_data))
    
    # 2. Final Results Section
    content_sections.append(generate_final_results_section(nipt_data))
    
    # 3. Trisomy Results Section
    content_sections.append(generate_trisomy_results_section(nipt_data))
    
    # 4. Trisomy Details Section
    content_sections.append(generate_trisomy_details_section(nipt_data, sample_id))
    
    # 5. Microdeletion Section
    content_sections.append(generate_microdeletion_section(nipt_data, sample_id))
    
    # 6. Quality Control Section
    content_sections.append(generate_quality_control_section(data, sample_id))
    
    # Combine all content
    content = '\n'.join(content_sections)
    
    # Generate timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Generate complete HTML
    template = generate_html_template()
    html_content = template.format(
        sample_id=sample_id,
        algorithm_version=algorithm_version,
        content=content,
        timestamp=timestamp,
        version=__version__
    )
    logger.info("HTML content generated successfully")
    
    # Save HTML file
    html_file = os.path.join(output_dir, f"{sample_id}_report.html")
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    logger.info(f"HTML report saved: {html_file}")
    logger.info(f"Expected image directory: {sample_dir}")
    
    return html_file

def generate_nipt_html_report(json_file_path, output_dir):
    """Main function to generate NIPT HTML report"""
    try:
        html_file = generate_html_report(json_file_path, output_dir)
        return html_file
    except Exception as e:
        logger.error(f"Error generating HTML report: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-json_file', dest='json_file', type=lambda x: file_check(parser, x), 
                       help='Input JSON file path')
    parser.add_argument('-output_dir', dest='output_dir', type=str,
                       help='Output directory for HTML file')
    parser.add_argument('-sample_name', dest='sample_name', 
                       help='Sample name (optional, will be extracted from JSON if not provided)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Generate HTML report
    try:
        html_file = generate_html_report(args.json_file, args.output_dir)
        logger.info(f"HTML report generated successfully: {html_file}")
        logger.info(f"Location: {os.path.abspath(html_file)}")
        
        # Optional: Open in browser
        try:
            import webbrowser
            file_url = f"file://{os.path.abspath(html_file)}"
            logger.info(f"🌐 Opening in browser: {file_url}")
            webbrowser.open(file_url)
        except:
            logger.error("💡 To view the report, open the HTML file in your web browser")
            
    except Exception as e:
        logger.error(f"❌ Error generating HTML report: {e}")
        return 1
    
    return 0

if __name__ == '__main__':
    exit(main())
