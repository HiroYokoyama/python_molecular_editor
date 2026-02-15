try:
    from PyQt6 import sip
except ImportError:
    import sip


def sip_isdeleted_safe(obj):
    """
    Safely check if a PyQt object has been deleted by C++.
    Returns True if deleted, False otherwise.
    """
    try:
        if obj is None:
            return True
        return sip.isdeleted(obj)
    except Exception:
        return False
