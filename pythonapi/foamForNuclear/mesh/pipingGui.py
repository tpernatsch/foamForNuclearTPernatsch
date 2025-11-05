
from PyQt5 import QtCore, QtGui, QtWidgets
import json

PORT_RADIUS = 6
HIT_RADIUS = 12
NODE_MIN_W_M = 0.05
NODE_MIN_H_M = 0.01
PX_PER_M = 100.0
PORT_OFFSET = 10.0

def clamp(v, lo, hi):
    return max(lo, min(hi, v))

def m_to_px(m): return m * PX_PER_M
def px_to_m(px): return px / PX_PER_M

class PortItem(QtWidgets.QGraphicsEllipseItem):
    def __init__(self, parent, side_id):
        r = PORT_RADIUS
        super().__init__(-r, -r, 2*r, 2*r, parent)
        self.setBrush(QtGui.QBrush(QtGui.QColor('white')))
        self.setPen(QtGui.QPen(QtGui.QColor('#111'), 1))
        self.setFlag(QtWidgets.QGraphicsItem.ItemIgnoresTransformations, True)
        self.side_id = side_id
        self.setZValue(3)
        self.setAcceptHoverEvents(True)
        self._hover = False
        self.edges = []

    def addEdge(self, edge):
        if edge not in self.edges:
            self.edges.append(edge)

    def removeEdge(self, edge):
        if edge in self.edges:
            self.edges.remove(edge)

    def hoverEnterEvent(self, e):
        self._hover = True
        self.update()
        super().hoverEnterEvent(e)

    def hoverLeaveEvent(self, e):
        self._hover = False
        self.update()
        super().hoverLeaveEvent(e)

    def paint(self, painter, option, widget):
        painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
        painter.setPen(QtGui.QPen(QtGui.QColor('#111'), 1))
        painter.setBrush(QtGui.QBrush(QtGui.QColor('#fff') if not self._hover else QtGui.QColor('#e0e7ff')))
        painter.drawEllipse(self.rect())

    def shape(self):
        r = PORT_RADIUS if self._hover else HIT_RADIUS
        path = QtGui.QPainterPath()
        path.addEllipse(QtCore.QRectF(-r, -r, 2*r, 2*r))
        return path

    def sceneAnchorPos(self):
        return self.mapToScene(self.rect().center())


class NodeItem(QtWidgets.QGraphicsRectItem):
    def __init__(self, nid, x=0, y=0, length_m=1.0, width_m=0.1, label='Pipe', color='#4f46e5', rotation=0, cells=1, kind='pipe'):
        self.kind = kind
        self.isPipe = self.kind == "pipe"
        super().__init__(0, 0, m_to_px(length_m), m_to_px(width_m) if self.isPipe else m_to_px(length_m))
        self.nid = nid
        self.length_m = float(max(length_m, NODE_MIN_W_M))
        self.width_m = float(max(width_m, NODE_MIN_H_M)) if self.isPipe else self.length_m
        self.cells = max(int(cells), 1)
        self.setPos(x, y)
        self.brushColor = QtGui.QColor(color)
        self.setPen(QtGui.QPen(QtGui.QColor('#1f2937'), 1))
        self.setFlags(QtWidgets.QGraphicsItem.ItemIsMovable | QtWidgets.QGraphicsItem.ItemIsSelectable | QtWidgets.QGraphicsItem.ItemIsFocusable | QtWidgets.QGraphicsItem.ItemSendsGeometryChanges)
        self.setTransformOriginPoint(self.rect().center())
        self.setRotation(rotation)
        self.labelItem = QtWidgets.QGraphicsTextItem(label, self)
        self.labelItem.setDefaultTextColor(QtGui.QColor('white'))
        self.labelItem.setTextInteractionFlags(QtCore.Qt.NoTextInteraction)
        self.labelItem.setZValue(2)
        # self.ports = {k: PortItem(self, k) for k in ('t','r','b','l')}
        self.ports = {k: PortItem(self, k) for k in ('r','l')}
        self._layout()
        self.updateToolTip()

    def updateToolTip(self):
        self.setToolTip(f"{self.labelItem.toPlainText()}\nType: {self.kind}\nLength: {self.length_m:.3f} m\nWidth: {self.width_m:.3f} m\nCells: {self.cells}")

    def _layout(self):
        r = self.rect()
        self.labelItem.setPos(r.center().x() - self.labelItem.boundingRect().width()/2, r.center().y() - self.labelItem.boundingRect().height()/2)
        self.ports['l'].setPos(r.left(), r.center().y())
        if (self.kind == "pump"):
            self.ports['r'].setPos(r.right(), r.center().y()-r.height()/4)
        elif (self.kind == "valve" or self.kind == "pipe"):
            self.ports['r'].setPos(r.right(), r.center().y())
            # self.ports['t'].setPos(r.center().x(), r.top())
            # self.ports['b'].setPos(r.center().x(), r.bottom())
        self.updateConnectedEdges()

    def updateConnectedEdges(self):
        for p in self.ports.values():
            for e in list(p.edges):
                e.updatePath()

    def paint(self, painter, option, widget):
        painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
        r = self.rect()
        painter.setPen(self.pen())
        painter.setBrush(self.brushColor)
        if self.kind == 'pipe':
            painter.drawRect(r)
        elif self.kind == 'valve':
            path = QtGui.QPainterPath()
            path.moveTo(r.left(), r.top())
            path.lineTo(r.left(), r.bottom())
            path.lineTo(r.right(), r.top())
            path.lineTo(r.right(), r.bottom())
            path.closeSubpath()
            painter.drawPath(path)
        elif self.kind == 'pump':
            path = QtGui.QPainterPath()
            path.moveTo(r.center().x(), r.top())
            path.lineTo(r.right(), r.top())
            path.lineTo(r.right(), r.center().y())
            path.arcTo(r.top(), r.left(), r.width(), r.height(), 0, -270)
            path.closeSubpath()
            painter.drawPath(path)
        if self.isSelected():
            pen = QtGui.QPen(QtGui.QColor('#2563eb'), 2, QtCore.Qt.DashLine)
            pen.setCosmetic(True)
            painter.setPen(pen)
            painter.setBrush(QtCore.Qt.NoBrush)
            painter.drawRect(r)

    def setColor(self, hexcolor):
        self.brushColor = QtGui.QColor(hexcolor)
        self.update()

    def setLabel(self, text):
        self.labelItem.setPlainText(text)
        self.updateToolTip()
        self._layout()

    def setCells(self, n):
        self.cells = max(int(n), 1)
        self.updateToolTip()

    def setLengthWidthMeters(self, length_m, width_m):
        self.length_m = max(float(length_m), NODE_MIN_W_M)
        self.width_m = max(float(width_m), NODE_MIN_H_M) if self.isPipe else self.length_m
        w_px = m_to_px(self.length_m)
        h_px = m_to_px(self.width_m)
        self.prepareGeometryChange()
        self.setRect(QtCore.QRectF(0, 0, w_px, h_px))
        self.setTransformOriginPoint(self.rect().center())
        self._layout()
        self.updateToolTip()

    def itemChange(self, change, value):
        if change in (QtWidgets.QGraphicsItem.ItemPositionChange, QtWidgets.QGraphicsItem.ItemTransformChange):
            self.updateConnectedEdges()
        out = super().itemChange(change, value)
        if change in (QtWidgets.QGraphicsItem.ItemPositionHasChanged, QtWidgets.QGraphicsItem.ItemTransformHasChanged):
            self.updateConnectedEdges()
        return out

    def toDict(self):
        r = self.rect()
        return {
            'id': self.nid,
            'type': self.kind,
            'position': {'x': float(self.pos().x()), 'y': float(self.pos().y())},
            'style': {'width_px': float(r.width()), 'height_px': float(r.height())},
            'data': {
                'label': self.labelItem.toPlainText(),
                'color': self.brushColor.name(),
                'rotation': float(self.rotation()),
                'length_m': float(self.length_m),
                'width_m': float(self.width_m),
                'cells': int(self.cells)
            }
        }


def _offset_from_side(pt, side):
    if side == 'l':
        return pt + QtCore.QPointF(-PORT_OFFSET, 0)
    if side == 'r':
        return pt + QtCore.QPointF(PORT_OFFSET, 0)
    if side == 't':
        return pt + QtCore.QPointF(0, -PORT_OFFSET)
    return pt + QtCore.QPointF(0, PORT_OFFSET)


class EdgeItem(QtWidgets.QGraphicsPathItem):
    def __init__(self, sourcePort, targetPort):
        super().__init__()
        self.sourcePort = sourcePort
        self.targetPort = targetPort
        self.setPen(QtGui.QPen(QtGui.QColor('#111'), 2))
        self.setZValue(1)
        self.setAcceptedMouseButtons(QtCore.Qt.LeftButton)
        self.setFlags(QtWidgets.QGraphicsItem.ItemIsSelectable)
        self.updatePath()
        self.sourcePort.addEdge(self)
        self.targetPort.addEdge(self)

    def updatePath(self):
        p1 = self.sourcePort.sceneAnchorPos()
        p2 = self.targetPort.sceneAnchorPos()
        # p1o = _offset_from_side(p1, self.sourcePort.side_id)
        # p2o = _offset_from_side(p2, self.targetPort.side_id)
        path = QtGui.QPainterPath(p1)
        # path.lineTo(p1o)
        # if abs(p1o.x() - p2o.x()) < 1e-3 or self.sourcePort.side_id in ('l','r') and self.targetPort.side_id in ('l','r'):
        #     mid = QtCore.QPointF(p2o.x(), p1o.y())
        # elif abs(p1o.y() - p2o.y()) < 1e-3 or self.sourcePort.side_id in ('t','b') and self.targetPort.side_id in ('t','b'):
        #     mid = QtCore.QPointF(p1o.x(), p2o.y())
        # else:
        #     mid = QtCore.QPointF(p2o.x(), p1o.y())
        if abs(p1.x() - p2.x()) < 1e-3 or self.sourcePort.side_id in ('l','r') and self.targetPort.side_id in ('l','r'):
            mid = QtCore.QPointF(p2.x(), p1.y())
        elif abs(p1.y() - p2.y()) < 1e-3 or self.sourcePort.side_id in ('t','b') and self.targetPort.side_id in ('t','b'):
            mid = QtCore.QPointF(p1.x(), p2.y())
        else:
            mid = QtCore.QPointF(p2.x(), p1.y())
        path.lineTo(mid)
        # path.lineTo(p2o)
        path.lineTo(p2)
        self.setPath(path)

    def cleanup(self):
        if hasattr(self, 'sourcePort') and self.sourcePort:
            self.sourcePort.removeEdge(self)
        if hasattr(self, 'targetPort') and self.targetPort:
            self.targetPort.removeEdge(self)
        self.sourcePort = None
        self.targetPort = None


class Scene(QtWidgets.QGraphicsScene):
    edgeCreated = QtCore.pyqtSignal(object)
    def __init__(self):
        super().__init__()
        self.tempLine = None
        self.dragSourcePort = None
        self.setSceneRect(0, 0, 4000, 3000)

    def _findPortAt(self, pos):
        for it in self.items(pos):
            if isinstance(it, PortItem):
                return it
            if isinstance(it, NodeItem):
                best_p = None
                best_d = 1e18
                for p in it.ports.values():
                    d = (p.sceneAnchorPos() - pos).manhattanLength()
                    if d < best_d:
                        best_d = d
                        best_p = p
                if best_p and best_d <= HIT_RADIUS:
                    return best_p
        return None

    def mousePressEvent(self, e):
        item = self._findPortAt(e.scenePos())
        if isinstance(item, PortItem) and e.button() == QtCore.Qt.LeftButton:
            self.dragSourcePort = item
            self.tempLine = QtWidgets.QGraphicsPathItem()
            self.tempLine.setPen(QtGui.QPen(QtGui.QColor('#555'), 1, QtCore.Qt.DashLine))
            self.tempLine.setZValue(0)
            self.tempLine.setAcceptedMouseButtons(QtCore.Qt.NoButton)
            self.addItem(self.tempLine)
            p = item.sceneAnchorPos()
            path = QtGui.QPainterPath(p)
            path.lineTo(p)
            self.tempLine.setPath(path)
            e.accept()
            return
        super().mousePressEvent(e)

    def mouseMoveEvent(self, e):
        if self.tempLine and self.dragSourcePort:
            p1 = self.dragSourcePort.sceneAnchorPos()
            p2 = e.scenePos()
            path = QtGui.QPainterPath(p1)
            path.lineTo(p2)
            self.tempLine.setPath(path)
            e.accept()
            return
        super().mouseMoveEvent(e)

    def mouseReleaseEvent(self, e):
        if self.tempLine and self.dragSourcePort:
            item = self._findPortAt(e.scenePos())
            if isinstance(item, PortItem) and item is not self.dragSourcePort:
                self.edgeCreated.emit((self.dragSourcePort, item))
            self.removeItem(self.tempLine)
            self.tempLine = None
            self.dragSourcePort = None
            e.accept()
            return
        super().mouseReleaseEvent(e)


class PropertyPanel(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._node = None
        self._block = False
        self.buildUI()

    def buildUI(self):
        form = QtWidgets.QFormLayout(self)
        self.nameEdit = QtWidgets.QLineEdit()
        self.colorBtn = QtWidgets.QPushButton("Pick…")
        self.rotSpin = QtWidgets.QDoubleSpinBox(); self.rotSpin.setRange(-180.0, 180.0); self.rotSpin.setSingleStep(1.0)
        self.lenSpin = QtWidgets.QDoubleSpinBox(); self.lenSpin.setRange(NODE_MIN_W_M, 1e9); self.lenSpin.setDecimals(3); self.lenSpin.setSuffix(" m")
        self.wSpin = QtWidgets.QDoubleSpinBox(); self.wSpin.setRange(NODE_MIN_H_M, 1e6); self.wSpin.setDecimals(3); self.wSpin.setSuffix(" m")
        self.cellsSpin = QtWidgets.QSpinBox(); self.cellsSpin.setRange(1, 10**9)
        form.addRow("Name", self.nameEdit)
        form.addRow("Color", self.colorBtn)
        form.addRow("Rotation", self.rotSpin)
        form.addRow("Pipe length", self.lenSpin)
        form.addRow("Pipe width", self.wSpin)
        form.addRow("Cells", self.cellsSpin)
        self.setDisabled(True)
        self.nameEdit.textEdited.connect(self.onName)
        self.colorBtn.clicked.connect(self.onColor)
        self.rotSpin.valueChanged.connect(self.onRot)
        self.lenSpin.valueChanged.connect(self.onLen)
        self.wSpin.valueChanged.connect(self.onW)
        self.cellsSpin.valueChanged.connect(self.onCells)

    def bindNode(self, node):
        self._node = node
        self._block = True
        try:
            if node is None:
                self.setDisabled(True)
                self.nameEdit.clear()
                self.rotSpin.setValue(0.0)
                self.lenSpin.setValue(NODE_MIN_W_M)
                self.wSpin.setValue(NODE_MIN_H_M)
                self.cellsSpin.setValue(1)
            else:
                self.setDisabled(False)
                self.nameEdit.setText(node.labelItem.toPlainText())
                self.rotSpin.setValue(node.rotation())
                self.lenSpin.setValue(node.length_m)
                self.wSpin.setValue(node.width_m)
                self.wSpin.setDisabled(not node.isPipe)
                self.cellsSpin.setValue(node.cells)
        finally:
            self._block = False

    def onName(self, s):
        if self._block or not self._node: return
        self._node.setLabel(s)
    def onColor(self):
        if self._block or not self._node: return
        col = QtWidgets.QColorDialog.getColor(parent=self)
        if col.isValid():
            self._node.setColor(col.name())
            self._node.updateConnectedEdges()
    def onRot(self, v):
        if self._block or not self._node: return
        self._node.setRotation(v)
        self._node.updateConnectedEdges()
    def onLen(self, v):
        if self._block or not self._node: return
        self._node.setLengthWidthMeters(v, self.wSpin.value())
    def onW(self, v):
        if self._block or not self._node: return
        self._node.setLengthWidthMeters(self.lenSpin.value(), v)
    def onCells(self, v):
        if self._block or not self._node: return
        self._node.setCells(v)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("foamForNuclear Plant Diagrammer")
        self.resize(1400, 850)
        self.scene = Scene()
        self.view = QtWidgets.QGraphicsView(self.scene)
        self.view.setRenderHints(QtGui.QPainter.Antialiasing | QtGui.QPainter.TextAntialiasing)
        self.view.setDragMode(QtWidgets.QGraphicsView.RubberBandDrag)
        self.setCentralWidget(self.view)
        self.toolbar = self.addToolBar("Tools")
        self._buildToolbar()
        self.nodes = {}
        self.edges = []
        self.scene.edgeCreated.connect(self._onEdgeCreated)
        self.scene.selectionChanged.connect(self._onSelectionChanged)
        self.propDock = QtWidgets.QDockWidget("Properties", self)
        self.propDock.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable | QtWidgets.QDockWidget.DockWidgetFloatable)
        self.propPanel = PropertyPanel(self)
        self.propDock.setWidget(self.propPanel)
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.propDock)

    def _buildToolbar(self):
        addPipe = QtWidgets.QAction("Pipe", self)
        addValve = QtWidgets.QAction("Valve", self)
        addPump = QtWidgets.QAction("Pump", self)
        imp = QtWidgets.QAction("Import", self)
        exp = QtWidgets.QAction("Export", self)
        clear = QtWidgets.QAction("New", self)
        color = QtWidgets.QAction("Color", self)
        rotate_l = QtWidgets.QAction("⟲", self)
        rotate_r = QtWidgets.QAction("⟳", self)
        delete_act = QtWidgets.QAction("Delete", self)
        delete_act.setShortcut(QtGui.QKeySequence.Delete)
        self.toolbar.addActions([addPipe, addValve, addPump])
        self.toolbar.addSeparator()
        self.toolbar.addActions([color, rotate_l, rotate_r, delete_act])
        self.toolbar.addSeparator()
        self.toolbar.addActions([imp, exp, clear])
        addPipe.triggered.connect(lambda: self.addNode('Pipe', '#0ea5e9', 'pipe'))
        addValve.triggered.connect(lambda: self.addNode('Valve', '#f97316', 'valve'))
        addPump.triggered.connect(lambda: self.addNode('Pump', '#22c55e', 'pump'))
        color.triggered.connect(self.changeColorSelected)
        rotate_l.triggered.connect(lambda: self.rotateSelected(-15))
        rotate_r.triggered.connect(lambda: self.rotateSelected(15))
        delete_act.triggered.connect(self.deleteSelected)
        imp.triggered.connect(self.importJSON)
        exp.triggered.connect(self.exportJSON)
        clear.triggered.connect(self.clearAll)

    def sceneCenter(self):
        return self.view.mapToScene(self.view.viewport().rect().center())

    def addNode(self, label, color, kind):
        nid = f"n{max([int(e.replace('n', '')) for e in self.nodes.keys()]+[0])+1}"
        center = self.sceneCenter()
        L = 1.0
        W = 0.1
        node = NodeItem(nid, center.x() - m_to_px(L)/2, center.y() - m_to_px(W)/2, L, W, label=label, color=color, rotation=0, cells=1, kind=kind)
        self.scene.addItem(node)
        self.nodes[nid] = node
        self.scene.clearSelection()
        node.setSelected(True)
        self._onSelectionChanged()

    def _onEdgeCreated(self, ports):
        s, t = ports
        edge = EdgeItem(s, t)
        self.scene.addItem(edge)
        self.edges.append(edge)

    def _onSelectionChanged(self):
        items = self.scene.selectedItems()
        node = next((i for i in items if isinstance(i, NodeItem)), None)
        self.propPanel.bindNode(node)

    def changeColorSelected(self):
        col = QtWidgets.QColorDialog.getColor(parent=self)
        if not col.isValid():
            return
        for item in self.scene.selectedItems():
            if isinstance(item, NodeItem):
                item.setColor(col.name())
                item.updateConnectedEdges()
        self._onSelectionChanged()

    def rotateSelected(self, delta):
        for item in self.scene.selectedItems():
            if isinstance(item, NodeItem):
                item.setRotation(item.rotation() + delta)
                item.updateConnectedEdges()
        self._onSelectionChanged()

    def deleteSelected(self):
        for item in list(self.scene.selectedItems()):
            if isinstance(item, EdgeItem):
                item.cleanup()
                if item in self.edges:
                    self.edges.remove(item)
                self.scene.removeItem(item)
            elif isinstance(item, NodeItem):
                for p in item.ports.values():
                    for e in list(p.edges):
                        e.cleanup()
                        if e in self.edges:
                            self.edges.remove(e)
                        self.scene.removeItem(e)
                nid_to_delete = [k for k, v in self.nodes.items() if v is item]
                for k in nid_to_delete:
                    del self.nodes[k]
                self.scene.removeItem(item)
        self._onSelectionChanged()

    def clearAll(self):
        self.scene.clear()
        self.nodes = {}
        self.edges = []
        self._onSelectionChanged()

    def exportJSON(self):
        data = {'version': 1, 'nodes': [], 'edges': []}
        for node in self.nodes.values():
            data['nodes'].append(node.toDict())
        for e in self.edges:
            def nodeIdForPort(p):
                parent = p.parentItem()
                return parent.nid if isinstance(parent, NodeItem) else None
            data['edges'].append({'source': nodeIdForPort(e.sourcePort), 'sourceHandle': e.sourcePort.side_id, 'target': nodeIdForPort(e.targetPort), 'targetHandle': e.targetPort.side_id, 'type': 'orthogonal'})
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export JSON", "diagram.json", "JSON (*.json)")
        if path:
            with open(path, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2)

    def importJSON(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Import JSON", "", "JSON (*.json)")
        if not path:
            return
        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        self.clearAll()
        id_to_node = {}
        for n in data.get('nodes', []):
            pos = n.get('position', {})
            d = n.get('data', {})
            L = d.get('length_m', px_to_m(n.get('style', {}).get('width_px', m_to_px(1.0))))
            W = d.get('width_m', px_to_m(n.get('style', {}).get('height_px', m_to_px(0.1))))
            C = int(d.get('cells', 1))
            kind = n.get('type', d.get('kind', 'pipe'))
            node = NodeItem(n.get('id', f"n{len(self.nodes)+1}"), pos.get('x', 0), pos.get('y', 0), L, W, label=d.get('label', 'Pipe'), color=d.get('color', '#4f46e5'), rotation=d.get('rotation', 0), cells=C, kind=kind)
            self.scene.addItem(node)
            self.nodes[node.nid] = node
            id_to_node[node.nid] = node
        for e in data.get('edges', []):
            s = id_to_node.get(e.get('source'))
            t = id_to_node.get(e.get('target'))
            if s and t:
                sp = s.ports.get(e.get('sourceHandle', 'r'), s.ports['r'])
                tp = t.ports.get(e.get('targetHandle', 'l'), t.ports['l'])
                edge = EdgeItem(sp, tp)
                self.scene.addItem(edge)
                self.edges.append(edge)

    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        if event.modifiers() & QtCore.Qt.ControlModifier:
            f = 1.15 if delta > 0 else 1/1.15
            self.view.scale(f, f)
            event.accept()
        else:
            super().wheelEvent(event)


class PipingGUI:
    def __init__(self):
        self.app = None
        self.win = None

    def run(self):
        import sys
        QtWidgets.QApplication.setDesktopSettingsAware(False)
        from PyQt5 import QtCore as _QtCore
        _QtCore.QCoreApplication.setAttribute(_QtCore.Qt.AA_UseSoftwareOpenGL, True)
        self.app = QtWidgets.QApplication(sys.argv)
        self.app.setStyle("Fusion")
        self.win = MainWindow()
        self.win.importJSON()
        self.win.show()
        sys.exit(self.app.exec_())
