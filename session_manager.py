import json
from typing import Dict, Any, Optional
from datetime import datetime


class SessionManager:
    """Save and load RNA-seq analysis sessions as JSON."""
    
    SERIALIZABLE_KEYS = [
        "comparisons", "sample_metadata", "conditions",
        "padj_threshold", "lfc_threshold", "active_comparison",
        "analysis_complete", "current_step",
    ]
    
    PLATFORM_VERSION = "1.1.0"
    
    @staticmethod
    def save_session(session_state: dict) -> str:
        """Serialize session state to JSON string.
        
        Only saves lightweight config/settings, NOT DataFrames or figures.
        """
        data = {
            "_meta": {
                "platform_version": SessionManager.PLATFORM_VERSION,
                "saved_at": datetime.now().isoformat(),
            }
        }
        for key in SessionManager.SERIALIZABLE_KEYS:
            if key in session_state:
                val = session_state[key]
                # Convert tuples to lists for JSON
                if isinstance(val, (list, tuple)):
                    data[key] = _make_serializable(val)
                else:
                    data[key] = val
        return json.dumps(data, indent=2, default=str)
    
    @staticmethod
    def load_session(json_str: str) -> dict:
        """Deserialize JSON back to session state dict."""
        data = json.loads(json_str)
        # Remove meta, return rest
        data.pop("_meta", None)
        return data
    
    @staticmethod
    def get_session_summary(session_state: dict) -> dict:
        """Get human-readable summary of current session."""
        summary = {
            "analysis_complete": session_state.get("analysis_complete", False),
            "num_comparisons": len(session_state.get("comparisons", [])),
            "conditions": session_state.get("conditions", []),
            "thresholds": {
                "padj": session_state.get("padj_threshold", 0.05),
                "lfc": session_state.get("lfc_threshold", 1.0),
            },
        }
        meta = session_state.get("sample_metadata", {})
        summary["num_samples"] = len(meta)
        return summary


def _make_serializable(obj):
    """Recursively convert tuples to lists for JSON serialization."""
    if isinstance(obj, tuple):
        return [_make_serializable(item) for item in obj]
    elif isinstance(obj, list):
        return [_make_serializable(item) for item in obj]
    elif isinstance(obj, dict):
        return {k: _make_serializable(v) for k, v in obj.items()}
    return obj
