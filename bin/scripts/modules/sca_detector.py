"""
SCA Detection Configuration Manager

SCA detection을 위한 설정 파일 관리 및 detection 수행 모듈
"""

import json
import os
import logging

logger = logging.getLogger(__name__)

class SCADetector:
    """
    SCA Detection 설정 파일 통합 관리 클래스
    """
    
    def __init__(self, data_dir=None, labcode=None):
        """
        Args:
            data_dir (str): DATA_DIR 경로
            labcode (str): 실험실 코드 (예: 'cordlife')
        """
        self.data_dir = data_dir
        self.labcode = labcode
        self.configs = {}
        self.current_config = None
        self.current_type = None
        
        # 설정 파일 경로 구성
        if data_dir and labcode:
            self.config_base_dir = os.path.join(data_dir, 'refs', labcode, 'EZD')
        else:
            self.config_base_dir = None
    
    def _get_config_file_path(self, config_type):
        """
        설정 타입에 따른 파일 경로 생성
        
        Args:
            config_type (str): 'orig', 'fetus', 'mom' 중 선택
        
        Returns:
            str: 설정 파일 전체 경로
        """
        if self.config_base_dir is None:
            # 기본 경로 (현재 디렉토리)
            config_files = {
                'orig': 'sca_config.json',
                'fetus': 'sca_config.json', 
                'mom': 'sca_config.json'
            }
            return config_files.get(config_type)
        
        # DATA_DIR/refs/<labcode>/EZD/<config_type>/sca_config.json
        config_dir = os.path.join(self.config_base_dir, config_type)
        config_file = os.path.join(config_dir, 'sca_config.json')
        
        return config_file
    
    @classmethod
    def create_default_configs(cls, data_dir, labcode):
        """
        기본 설정 파일들을 적절한 디렉토리에 생성
        
        Args:
            data_dir (str): DATA_DIR 경로
            labcode (str): 실험실 코드
        """
        base_dir = os.path.join(data_dir, 'refs', labcode, 'EZD')
        
        # 설정 데이터들
        configs = {
            'orig': {
                "config_info": {
                    "name": "Original Reference SCA Detection",
                    "type": "orig",
                    "description": "Original population reference data for SCA detection",
                    "version": "1.0",
                    "created_date": "2025-05-31"
                },
                "sca_detection": {
                    "male": {
                        "slope": -0.08367685038621694,
                        "intercept": 0.4876543737788628,
                        "ur_x_threshold": 5.2,
                        "margin": 0.005
                    },
                    "female": {
                        "xo_z_threshold": -6.0,
                        "xxx_z_threshold": 4.5,
                        "ur_x_low": 5.35,
                        "ur_x_high": 5.45,
                        "z_normal_low": -3.0,
                        "z_normal_high": 1.0,
                        "xo_ur_x_min": 4.9,
                        "xo_ur_x_max": 5.2,
                        "xxx_ur_x_min": 5.6,
                        "xxx_ur_x_max": 6.0
                    }
                }
            },
            'fetus': {
                "config_info": {
                    "name": "Fetal Reference SCA Detection",
                    "type": "fetus",
                    "description": "Fetal population reference data for SCA detection",
                    "version": "1.0",
                    "created_date": "2025-05-31"
                },
                "sca_detection": {
                    "male": {
                        "slope": -0.09061059750410849,
                        "intercept": 0.5229228546732532,
                        "ur_x_threshold": 5.2,
                        "margin": 0.005
                    },
                    "female": {
                        "xo_z_threshold": -6.5,
                        "xxx_z_threshold": 4.0,
                        "ur_x_low": 5.3,
                        "ur_x_high": 5.5,
                        "z_normal_low": -3.2,
                        "z_normal_high": 1.2,
                        "xo_ur_x_min": 4.8,
                        "xo_ur_x_max": 5.3,
                        "xxx_ur_x_min": 5.5,
                        "xxx_ur_x_max": 6.2
                    }
                }
            },
            'mom': {
                "config_info": {
                    "name": "Maternal Reference SCA Detection",
                    "type": "mom",
                    "description": "Maternal population reference data for SCA detection",
                    "version": "1.0",
                    "created_date": "2025-05-31"
                },
                "sca_detection": {
                    "male": {
                        "enabled": False,
                        "reason": "Mom reference does not require male SCA detection"
                    },
                    "female": {
                        "xo_z_threshold": -5.5,
                        "xxx_z_threshold": 5.0,
                        "ur_x_low": 5.4,
                        "ur_x_high": 5.6,
                        "z_normal_low": -2.5,
                        "z_normal_high": 1.5,
                        "xo_ur_x_min": 5.0,
                        "xo_ur_x_max": 5.3,
                        "xxx_ur_x_min": 5.7,
                        "xxx_ur_x_max": 6.5
                    }
                }
            }
        }
        
        # 디렉토리 생성 및 파일 저장
        for config_type, config_data in configs.items():
            config_dir = os.path.join(base_dir, config_type)
            os.makedirs(config_dir, exist_ok=True)
            
            config_file = os.path.join(config_dir, 'sca_config.json')
            
            with open(config_file, 'w', encoding='utf-8') as f:
                json.dump(config_data, f, indent=4, ensure_ascii=False)
            
            logger.info(f"Config file generation: {config_file}")
        
        logger.info(f"All SCA config files were generated in {base_dir}.")
    
    def load_config(self, config_type='orig'):
        """
        특정 타입의 설정 파일 로드
        
        Args:
            config_type (str): 'orig', 'fetus', 'mom' 중 선택
        
        Returns:
            dict: 로드된 설정 데이터
        """
        if config_type not in ['orig', 'fetus', 'mom']:
            logger.info(f"Not supported config type: {config_type}")
            return None
        
        config_file = self._get_config_file_path(config_type)
        
        if config_file is None:
            logger.info(f"Cannot find Config file : {config_type}")
            return None
        
        try:
            logger.info(f"{config_file}")
            with open(config_file, 'r', encoding='utf-8') as f:
                config = json.load(f)
            
            self.configs[config_type] = config
            self.current_config = config
            self.current_type = config_type
            
            config_name = config.get('config_info', {}).get('name', f'{config_type} config')
            logger.info(f"SCA confi loading completed!: {config_name}")
            return config
            
        except FileNotFoundError:
            logger.info(f"Config file not found : {config_file}")
            logger.info(f"SCADetector.create_default_configs('{self.data_dir}', '{self.labcode}')")
            return None

        except json.JSONDecodeError as e:
            logger.info(f"JSON format error : {e}")
            return None
    
    def load_all_configs(self):
        """모든 설정 파일 로드"""
        config_types = ['orig', 'fetus', 'mom']
        loaded_count = 0
        
        for config_type in config_types:
            if self.load_config(config_type):
                loaded_count += 1
        
        logger.info(f"{loaded_count} SCA config files were loaded")
        return loaded_count
    
    def get_male_params(self, config_type='orig'):
        """
        Male SCA detection 파라미터 추출
        
        Args:
            config_type (str): 'orig', 'fetus', 'mom' 중 선택
        
        Returns:
            dict: male detection 파라미터들
        """
        if config_type not in self.configs:
            self.load_config(config_type)
        
        if config_type not in self.configs:
            return None
        
        config = self.configs[config_type]
        male_config = config['sca_detection']['male']
        
        # mom 설정은 male detection이 비활성화됨
        if not male_config.get('enabled', True):
            logger.info(f"male SCA detection was not supported in {config_type}")
            return None
        
        return {
            'slope': male_config['slope'],
            'intercept': male_config['intercept'],
            'ur_x_threshold': male_config['ur_x_threshold'],
            'margin': male_config['margin'],
            'description': male_config['description'],
            'sample_count': male_config.get('sample_count', 0)
        }
    
    def get_female_params(self, config_type='orig'):
        """
        Female SCA detection 파라미터 추출
        
        Args:
            config_type (str): 'orig', 'fetus', 'mom' 중 선택
        
        Returns:
            dict: female detection 파라미터들
        """
        if config_type not in self.configs:
            self.load_config(config_type)
        
        if config_type not in self.configs:
            return None
        
        config = self.configs[config_type]
        female_config = config['sca_detection']['female']
        
        return {
            'xo_z_threshold': female_config['xo_z_threshold'],
            'xxx_z_threshold': female_config['xxx_z_threshold'],
            'ur_x_low': female_config['ur_x_low'],
            'ur_x_high': female_config['ur_x_high'],
            'z_normal_low': female_config['z_normal_low'],
            'z_normal_high': female_config['z_normal_high'],
            'xo_ur_x_min': female_config['xo_ur_x_min'],
            'xo_ur_x_max': female_config['xo_ur_x_max'],
            'xxx_ur_x_min': female_config['xxx_ur_x_min'],
            'xxx_ur_x_max': female_config['xxx_ur_x_max'],
            'description': female_config['description'],
            'sample_count': female_config.get('sample_count', 0)
        }
    
    def detect_male_sca(self, ur_x, ur_y, config_type='orig'):
        """
        Male SCA detection 수행
        
        Args:
            ur_x (float): UAR[X] 값
            ur_y (float): UAR[Y] 값 
            config_type (str): 'orig', 'fetus' 중 선택 ('mom'은 지원 안함)
        
        Returns:
            str: detection 결과
        """
        params = self.get_male_params(config_type)
        if params is None:
            return "설정 오류 또는 지원하지 않는 타입"
        
        # 경계선의 y값 계산
        boundary_y = params['slope'] * ur_x + params['intercept']
        
        # SCA 판정
        if ur_y > boundary_y + params['margin']:
            if ur_x <= params['ur_x_threshold']:
                return "XYY Detected"
            else:
                return "XXY Detected"
        elif ur_y > boundary_y:
            if ur_x <= params['ur_x_threshold']:
                return "XYY Suspected"
            else:
                return "XXY Suspected"
        else:
            return "Not Detected"

    def detect_female_sca(self, ur_x, z_score, config_type='orig'):

        params = self.get_female_params(config_type)
        if params is None:
            return "설정 오류"
        
        # UR_X 영역 분류
        def classify_ur_x_region(ur_x_val):
            if ur_x_val < params['xo_ur_x_max']:
                return 'XO_DETECTED'
            elif ur_x_val < params['ur_x_low']:
                return 'XO_SUSPECTED'
            elif ur_x_val <= params['ur_x_high']:
                return 'NORMAL'
            elif ur_x_val < params['xxx_ur_x_min']:
                return 'XXX_SUSPECTED'
            else:
                return 'XXX_DETECTED'
        
        # Z-score 영역 분류
        def classify_z_region(z_val):
            if z_val < params['xo_z_threshold']:
                return 'XO'
            elif z_val <= params['z_normal_high']:
                return 'NORMAL'
            elif z_val > params['xxx_z_threshold']:
                return 'XXX'
            else:
                return 'BORDERLINE'
        
        ur_x_region = classify_ur_x_region(ur_x)
        z_region = classify_z_region(z_score)
        
        # 결과 매핑 테이블
        result_map = {
            # (Z_region, UR_X_region): Result
            ('NORMAL', 'NORMAL'): "Not Detected",
            
            # XO cases
            ('XO', 'XO_DETECTED'): "XO Detected",
            ('XO', 'XO_SUSPECTED'): "XO Suspected",
            ('XO', 'NORMAL'): "XO Suspected",
            ('NORMAL', 'XO_DETECTED'): "XO Suspected",
            ('NORMAL', 'XO_SUSPECTED'): "XO Suspected",
            
            # XXX cases  
            ('XXX', 'XXX_DETECTED'): "XXX Detected",
            ('XXX', 'XXX_SUSPECTED'): "XXX Suspected",
            ('XXX', 'NORMAL'): "XXX Suspected",
            ('NORMAL', 'XXX_DETECTED'): "XXX Suspected",
            ('NORMAL', 'XXX_SUSPECTED'): "XXX Suspected",
            
            # Borderline cases
            ('BORDERLINE', 'NORMAL'): "Not Detected",
            ('BORDERLINE', 'XO_SUSPECTED'): "XO Suspected",
            ('BORDERLINE', 'XO_DETECTED'): "XO Suspected",
            ('BORDERLINE', 'XXX_SUSPECTED'): "XXX Suspected",
            ('BORDERLINE', 'XXX_DETECTED'): "XXX Suspected",
            
            # Cross cases (XO z-score + XXX ur_x or vice versa)
            ('XO', 'XXX_SUSPECTED'): "Not Detected",  # 상충되는 신호
            ('XO', 'XXX_DETECTED'): "Not Detected",   # 상충되는 신호
            ('XXX', 'XO_SUSPECTED'): "Not Detected",  # 상충되는 신호
            ('XXX', 'XO_DETECTED'): "Not Detected",   # 상충되는 신호
        }
        
        # 결과 반환
        key = (z_region, ur_x_region)
        return result_map.get(key, "Not Detected")

    def detect_female_sca_old(self, ur_x, z_score, config_type='orig'):
        """
        Female SCA detection 수행
        
        Args:
            ur_x (float): UAR[X] 값
            z_score (float): Z-score 값
            config_type (str): 'orig', 'fetus', 'mom' 중 선택
        
        Returns:
            str: detection 결과
        """
        params = self.get_female_params(config_type)
        if params is None:
            return "설정 오류"
        
        # XO Detection
        # xo_ur_x_min may drop much, so, don't consider this
        if (z_score < params['xo_z_threshold'] and 
            #params['xo_ur_x_min'] <= ur_x <= params['xo_ur_x_max']):
            ur_x <= params['xo_ur_x_max']):
            return "XO Detected"
        
        # XXX Detection  
        # xxx_ur_x_max may increase much, so, don't consider this
        if (z_score > params['xxx_z_threshold'] and
            #params['xxx_ur_x_min'] <= ur_x <= params['xxx_ur_x_max']):
            ur_x >= params['xxx_ur_x_min']):
            return "XXX Detected"
        
        # Normal
        if (params['z_normal_low'] <= z_score <= params['z_normal_high'] and
            params['ur_x_low'] <= ur_x <= params['ur_x_high']):
            return "Not Detected"
        
        return "Not Detected"
    
    def compare_configs(self, ur_x, ur_y=None, z_score=None, gender='male'):
        """
        여러 설정으로 detection 결과 비교
        
        Args:
            ur_x (float): UAR[X] 값
            ur_y (float): UAR[Y] 값 (male용)
            z_score (float): Z-score 값 (female용)
            gender (str): 'male' 또는 'female'
        
        Returns:
            dict: 설정별 결과
        """
        results = {}
        
        config_types = ['orig', 'fetus', 'mom'] if gender == 'female' else ['orig', 'fetus']
        
        for config_type in config_types:
            try:
                if gender == 'male' and ur_y is not None:
                    result = self.detect_male_sca(ur_x, ur_y, config_type)
                    params = self.get_male_params(config_type)
                    boundary_y = params['slope'] * ur_x + params['intercept'] if params else None
                    
                    results[config_type] = {
                        'result': result,
                        'boundary_y': boundary_y,
                        'distance': ur_y - boundary_y if boundary_y else None,
                        'config_name': self.configs.get(config_type, {}).get('config_info', {}).get('name', config_type)
                    }
                
                elif gender == 'female' and z_score is not None:
                    result = self.detect_female_sca(ur_x, z_score, config_type)
                    
                    results[config_type] = {
                        'result': result,
                        'ur_x': ur_x,
                        'z_score': z_score,
                        'config_name': self.configs.get(config_type, {}).get('config_info', {}).get('name', config_type)
                    }
                    
            except Exception as e:
                results[config_type] = {'result': f'오류: {e}'}
        
        return results
    
    def print_male_results(self, male_results, ur_x, ur_y):
        """Male SCA detection 결과를 깔끔하게 출력"""
        logger.info("=" * 60)
        logger.info(f"Male SCA Detection result")
        logger.info("=" * 60)
        logger.info(f"테스트 샘플: UAR[X] = {ur_x:.4f}, UAR[Y] = {ur_y:.4f}")
        logger.info("-" * 60)
        
        for config_type, result_info in male_results.items():
            result = result_info.get('result', 'N/A')
            boundary_y = result_info.get('boundary_y', None)
            distance = result_info.get('distance', None)
            config_name = result_info.get('config_name', config_type)
            
            logger.info(f"config: {config_type.upper():6} ({config_name})")
            logger.info(f"  result: {result}")
            
            if boundary_y is not None:
                logger.info(f"  boundary Y: {boundary_y:.6f}")
                
            if distance is not None:
                status = "UP" if distance > 0 else "DOWN"
                logger.info(f"  Distance from the boundary: {distance:+.6f} ({status})")
                

    def print_female_results(self, female_results, ur_x, z_score):
        """Female SCA detection 결과를 깔끔하게 출력"""
        logger.info("=" * 60)
        logger.info(f"Female SCA Detection result")
        logger.info("=" * 60)
        logger.info(f"테스트 샘플: UAR[X] = {ur_x:.4f}, Z-score = {z_score:.4f}")
        logger.info("-" * 60)
        
        for config_type, result_info in female_results.items():
            result = result_info.get('result', 'N/A')
            config_name = result_info.get('config_name', config_type)
            
            logger.info(f"설정: {config_type.upper():6} ({config_name})")
            logger.info(f"  결과: {result}")

    def analyze_config_differences(self, male_results):
        """설정별 결과 차이점 분석"""
        logger.info("=" * 60)
        logger.info("Config differences")
        logger.info("=" * 60)
        
        configs = list(male_results.keys())
        results = [male_results[config]['result'] for config in configs]
        
        # 결과가 모두 같은지 확인
        if len(set(results)) == 1:
            logger.info("✅ Same result among configs")
            logger.info(f"   Common results: {results[0]}")
        else:
            logger.info("⚠️  Different result among configs")
            for config in configs:
                result = male_results[config]['result']
                boundary_y = male_results[config].get('boundary_y', 'N/A')
                logger.info(f"   {config.upper()}: {result} (boundary: {boundary_y})")
        
        # 경계선 차이 분석
        boundary_values = []
        for config in configs:
            boundary_y = male_results[config].get('boundary_y')
            if boundary_y is not None:
                boundary_values.append((config, boundary_y))
        
        if len(boundary_values) >= 2:
            logger.info(f"\nGap value of boundary:")
            for i, (config1, boundary1) in enumerate(boundary_values):
                for config2, boundary2 in boundary_values[i+1:]:
                    diff = abs(boundary1 - boundary2)
                    logger.info(f"   {config1.upper()} vs {config2.upper()}: {diff:.6f} difference")
        

    def quick_male_analysis(self, ur_x, ur_y, config_type='orig'):
        """Male SCA 빠른 분석"""
        if not self.configs:
            self.load_all_configs()
        
        if config_type == 'all':
            male_results = self.compare_configs(ur_x, ur_y=ur_y, gender='male')
            self.logger.info_male_results(male_results, ur_x, ur_y)
            self.analyze_config_differences(male_results)
            return male_results
        else:
            result = self.detect_male_sca(ur_x, ur_y, config_type)
            logger.info(f"Male SCA ({config_type}): {result}")
            return result

    def quick_female_analysis(self, ur_x, z_score, config_type='orig'):
        """Female SCA 빠른 분석"""
        if not self.configs:
            self.load_all_configs()
        
        if config_type == 'all':
            female_results = self.compare_configs(ur_x, z_score=z_score, gender='female')
            self.logger.info_female_results(female_results, ur_x, z_score)
            return female_results
        else:
            result = self.detect_female_sca(ur_x, z_score, config_type)
            logger.info(f"Female SCA ({config_type}): {result}")
            return result

    def get_plot_colors(self, config_type='orig'):
        """plotting용 색상 설정 추출"""
        if config_type not in self.configs:
            self.load_config(config_type)
        
        if config_type not in self.configs:
            return self._default_colors()
        
        try:
            return self.configs[config_type]['plot_settings']['colors']
        except KeyError:
            return self._default_colors()
    
    def _default_colors(self):
        """기본 색상 설정"""
        return {
            "xy_normal": "blue",
            "xxy": "orange",
            "xyy": "green", 
            "xx_normal": "gray",
            "xo": "coral",
            "xxx": "gold",
            "boundary_line": "purple",
            "threshold_line": "red",
            "test_sample": "red"
        }

    def get_available_configs(self):
        """사용 가능한 설정 타입들 반환"""
        return list(self.configs.keys()) if self.configs else []

# 편의 함수들
def quick_detect_male(ur_x, ur_y, config_type='orig', config_dir='configs'):
    """
    빠른 Male SCA detection
    
    Args:
        ur_x (float): UAR[X] 값
        ur_y (float): UAR[Y] 값
        config_type (str): 'orig', 'fetus' 중 선택
        config_dir (str): 설정 파일 디렉토리
    
    Returns:
        str: detection 결과
    """
    manager = SCADetector(config_dir)
    return manager.quick_male_analysis(ur_x, ur_y, config_type)

# 테스트 코드
if __name__ == "__main__":
    logger.info("=== SCA Config Manager 테스트 ===")
    
    # 매니저 생성 및 테스트
    manager = SCADetector()
    manager.load_all_configs()
    
    # 빠른 테스트
    logger.info("\n1. 빠른 Male 테스트:")
    manager.quick_male_analysis(5.27, 0.033, 'all')
    
    logger.info("\n2. 빠른 Female 테스트:")
    manager.quick_female_analysis(5.4, -4.5, 'all')
    
    logger.info("\n3. 편의 함수 테스트:")
    result = quick_detect_male(5.27, 0.033, 'orig')
    logger.info(f"편의 함수 결과: {result}")
